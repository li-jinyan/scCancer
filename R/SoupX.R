#' Remove background contamination from count matrix
#'
#' After the level of background contamination has been estimated or specified for a channel, calculate the resulting corrected count matrix with background contamination removed.
#'
#' This essentially subtracts off the mean expected background counts for each gene, then redistributes any "unused" counts.  A count is unused if its subtraction has no effect.  For example, subtracting a count from a gene that has zero counts to begin with.
#'
#' As expression data is highly sparse at the single cell level, it is highly recommended that clustering information be provided to allow the subtraction method to share information between cells.  Without grouping cells into clusters, it is difficult (and usually impossible) to tell the difference between a count of 1 due to background contamination and a count of 1 due to endogenous expression.  This ambiguity is removed at the cluster level where counts can be aggregated across cells.  This information can then be propagated back to the individual cell level to provide a more accurate removal of contaminating counts.
#'
#' To provide clustering information, either set clustering on the SoupChannel object with \code{\link{setClusters}} or explicitly passing the \code{clusters} parameter.
#'
#' If \code{roundToInt=TRUE}, this function will round the result to integers.  That is, it will take the floor of the connected value and then round back up with probability equal to the fractional part of the number.
#'
#' The \code{method} parameter controls how the removal of counts in performed.  This should almost always be left at the default ('subtraction'), which iteratively subtracts counts from all genes as described above.  The 'soupOnly' method will use a p-value based estimation procedure to identify those genes that can be confidently identified as having endogenous expression and removes everything else (described in greater detail below).  Because this method either removes all or none of the expression for a gene in a cell, the correction procedure is much faster.  Finally, the 'multinomial' method explicitly maximises the multinomial likelihood for each cell.  This method gives essentially identical results as 'subtraction' and is considerably slower.
#'
#' In greater detail, the 'soupOnly' method is done by sorting genes within each cell by their p-value under the null of the expected soup fraction using a Poisson model.  So that genes that definitely do have a endogenous contribution are at the end of the list with p=0.  Those genes for which there is poor evidence of endogenous cell expression are removed, until we have removed approximately nUMIs*rho molecules.  The cut-off to prevent removal of genes above nUMIs*rho in each cell is achieved by calculating a separate p-value for the total number of counts removed to exceed nUMIs*rho, again using a Poisson model.  The two p-values are combined using Fisher's method and the cut-off is applied to the resulting combined p-value calculated using a chi-squared distribution with 4 degrees of freedom.
#'
#' @param sc A SoupChannel object.
#' @param clusters A vector of cluster IDs, named by cellIDs.  If NULL clusters auto-loaded from \code{sc}.  If FALSE, no clusters are used.  See details.
#' @param method Method to use for correction.  See details.  One of 'multinomial', 'soupOnly', or 'subtraction'
#' @param roundToInt Should the resulting matrix be rounded to integers?
#' @param verbose Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.
#' @param tol Allowed deviation from expected number of soup counts.  Don't change this.
#' @param pCut The p-value cut-off used when \code{method='soupOnly'}.
#' @param ... Passed to expandClusters.
#' @return A modified version of the table of counts, with background contamination removed.
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom stats rbinom pchisq
adjustCounts = function(sc,clusters=NULL,method=c('subtraction','soupOnly','multinomial'),roundToInt=FALSE,verbose=1,tol=1e-3,pCut=0.01,...){
  #####################
  # Parameter checking
  method = match.arg(method)
  if(!is(sc,'SoupChannel'))
    stop("sc must be an object of type SoupChannel")
  if(!'rho' %in% colnames(sc$metaData))
    stop("Contamination fractions must have already been calculated/set.")
  #Check the clusters parameter. If it's null, try and auto-fetch
  if(is.null(clusters)){
    if('clusters' %in% colnames(sc$metaData)){
      clusters = setNames(as.character(sc$metaData$clusters),rownames(sc$metaData))
    }else{
      warning("Clustering data not found.  Adjusting counts at cell level.  You will almost certainly get better results if you cluster data first.")
      clusters=FALSE
    }
  }
  ################################################
  # Recursive application for when using clusters
  if(clusters[1]!=FALSE){
    #Check we have coverage of everything
    if(!all(colnames(sc$toc) %in% names(clusters)))
      stop("Invalid cluster specification.  clusters must be a named vector with all column names in the table of counts appearing.")
    #OK proceed
    s = split(colnames(sc$toc),clusters[colnames(sc$toc)])
    tmp = sc
    #Create by cluster table of counts.  Must be sparse matrix
    tmp$toc = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
    tmp$toc = Matrix(tmp$toc,sparse=TRUE)
    tmp$metaData = data.frame(nUMIs = sapply(s,function(e) sum(sc$metaData[e,'nUMIs'])),
                              rho = sapply(s,function(e) sum(sc$metaData[e,'rho']*sc$metaData[e,'nUMIs'])/sum(sc$metaData[e,'nUMIs'])))
    #Recursively apply.  Could be done more neatly I'm sure.  Don't round to integer until we've re-expanded.
    out = adjustCounts(tmp,clusters=FALSE,method=method,roundToInt=FALSE,verbose=verbose,tol=tol,pCut=pCut)
    #This gives the corrected table, we want the soup counts
    out = tmp$toc - out
    #Finally re-expand the results back to single cell level
    out = expandClusters(out,sc$toc,clusters,sc$metaData$nUMIs*sc$metaData$rho,verbose=verbose,...)
    #And convert back to a corrected table of counts
    out = sc$toc - out
  }else{
    ##############################
    # Actual adjustment of counts
    if(method=='multinomial'){
      #Quick initial guess at best fit
      if(verbose>1)
        message("Initialising with subtraction method.")
      fitInit = sc$toc - adjustCounts(sc,clusters=FALSE,method='subtraction',roundToInt=TRUE)
      ps = sc$soupProfile$est
      out = list()
      #Loop over cells
      if(verbose>0){
        message(sprintf("Fitting multinomial distribution to %d cells/clusters.",ncol(sc$toc)))
        pb=initProgBar(1,ncol(sc$toc))
      }
      for(i in seq(ncol(sc$toc))){
        if(verbose>0)
          setTxtProgressBar(pb,i)
        #How many soup molecules do we expect for this cell?
        tgtN = round(sc$metaData$rho[i]*sc$metaData$nUMIs[i])
        #And what are the observational limits on which genes they can be assigned to
        lims = sc$toc[,i]
        #Initialise
        fit = fitInit[,i]
        while(TRUE){
          #Work out which we can increase
          increasable = fit<lims
          decreasable = fit>0
          #And what the likelihood gain/cost for changing them is
          delInc = log(ps[increasable]) - log(fit[increasable]+1)
          delDec = -log(ps[decreasable]) +log(fit[decreasable])
          #Decide which swap(s) would lead to best likelihood gain
          wInc = wIncAll = which(increasable)[which(delInc==max(delInc))]
          wDec = wDecAll = which(decreasable)[which(delDec==max(delDec))]
          nInc = length(wIncAll)
          nDec = length(wDecAll)
          if(nInc>1)
            wInc = sample(wIncAll,1)
          if(nDec>1)
            wDec = sample(wDecAll,1)
          #How many do we have
          curN = sum(fit)
          if(curN<tgtN){
            if(verbose>2)
              message(sprintf("# choices: nInc=%d nDec=%d, Under-allocated (%d of %d), increasing...",nInc,nDec,curN,tgtN))
            fit[wInc] = fit[wInc]+1
          }else if(curN>tgtN){
            if(verbose>2)
              message(sprintf("# choices: nInc=%d nDec=%d, Over-allocated (%d of %d), decreasing...",nInc,nDec,curN,tgtN))
            fit[wDec] = fit[wDec]-1
          }else{
            #Check if the swap will increase likelihood
            delTot = max(delInc)+max(delDec)
            #Three possibilites
            #1. del>0, keep going, we're improving the fit
            #2. del<0, stop, we won't get any better
            #3. del==0, accumulate all the ones that can go up/down.  As del==0 this is reversable, so want to share out "excess" counts between this group.
            if(verbose>2)
              message(sprintf("# choices: nInc=%d nDec=%d, Total log likelihood difference %s",nInc,nDec,delTot))
            if(delTot==0){
              #As the difference is zero, all movements between wInc and wDec are reversible.  So want to distribute evenly the "available" counts between those in the ambiguous set.
              #Take them away from those that presently have them
              fit[wDecAll] = fit[wDecAll]-1
              #And share them equally amongst those in bucket
              zeroBucket = unique(c(wIncAll,wDecAll))
              fit[zeroBucket] = fit[zeroBucket] + length(wDecAll)/length(zeroBucket)
              if(verbose>2)
                message(sprintf("Ambiguous final configuration. Shared %d reads between %d equally likely options",length(wDecAll),length(zeroBucket)))
              break
            }else if(delTot<0){
              #Minimum reached.
              if(verbose>2)
                message("Unique final configuration.")
              break
            }else{
              #Improvements to be made.
              fit[wInc] = fit[wInc]+1
              fit[wDec] = fit[wDec]-1
            }
          }
        }
        out[[i]]=fit
      }
      if(verbose>0)
        close(pb)
      out = do.call(cbind,out)
      out = as(out,'dgTMatrix')
      rownames(out) = rownames(sc$toc)
      colnames(out) = colnames(sc$toc)
      out = sc$toc - out
    }else if(method=='soupOnly'){
      if(verbose>0)
        message("Identifying and removing genes likely to be pure contamination in each cell.")
      #Start by calculating the p-value against the null of soup.
      out = as(sc$toc,'dgTMatrix')
      if(verbose>1)
        message("Calculating probability of each gene being soup")
      p = ppois(out@x-1,sc$metaData$nUMIs[out@j+1]*sc$soupProfile$est[out@i+1]*sc$metaData$rho[out@j+1],lower.tail=FALSE)
      #Order them by cell, then by p-value
      o = order(-(out@j+1),p,decreasing=TRUE)
      #Get the running total for removal.  Could probably make this faster with some tricks.
      if(verbose>1)
        message("Calculating probability of the next count being soup")
      s = split(o,out@j[o]+1)
      rTot = unlist(lapply(s,function(e) cumsum(out@x[e])),use.names=FALSE)
      #Now we need to get the soup probability vector as well.
      pSoup = ppois(rTot-out@x[o]-1,sc$metaData$nUMIs[out@j[o]+1]*sc$metaData$rho[out@j[o]+1],lower.tail=FALSE)
      if(verbose>1)
        message("Filtering table of counts")
      #Now construct the final probability vector
      pp = p[o]*pSoup
      q = pchisq(-2*log(pp),4,lower.tail=FALSE)
      #And we then want to drop anything meeting our termination criteria
      w = which(q<pCut)
      #Keep a list of dropped genes
      dropped = data.frame(cell=colnames(out)[out@j[o[-w]]+1],
                           gene=rownames(out)[out@i[o[-w]]+1],
                           cnt = out@x[o[-w]])
      if(verbose>2){
        message(sprintf("Most removed genes are:"))
        x = sort(table(dropped$gene)/ncol(out),decreasing=TRUE)
        print(x[seq_len(min(length(x),100))])
      }
      #Construct the corrected count matrix
      out = sparseMatrix(i=out@i[o[w]]+1,
                         j=out@j[o[w]]+1,
                         x=out@x[o[w]],
                         dims=dim(out),
                         dimnames=dimnames(out),
                         giveCsparse=FALSE)
    }else if(method=='subtraction'){
      #Create the final thing efficiently without making a big matrix
      out = as(sc$toc,'dgTMatrix')
      expSoupCnts = sc$metaData$nUMIs * sc$metaData$rho
      soupFrac = sc$soupProfile$est
      #Distribute counts according to the soup profile.  Could be made faster by not considering zeros, but eh.
      out = out - do.call(cbind,lapply(seq(ncol(out)),function(e) alloc(expSoupCnts[e],out[,e],soupFrac)))
      out = as(out,'dgTMatrix')
      #Iteratively remove until we've removed the expected number of counts
      #How many counts should we have left when we're done?
      #tgts = sc$metaData$nUMIs - expSoupCnts
      #Which genes do we still need to bother trying to remove counts from in which cells
      #toAdjust = seq_along(out@i)
      #if(verbose>0)
      #  message("Subtracting contaminating counts")
      #while(TRUE){
      #  #How many left to do?
      #  toDo = colSums(out)-tgts
      #  #Get the soup frac correction factors.
      #  tmp = rep(1,length(toDo))
      #  soupFracSum = sapply(split(soupFrac[out@i[toAdjust]+1],out@j[toAdjust]+1),sum)
      #  tmp[as.numeric(names(soupFracSum))]=soupFracSum
      #  toDo = toDo/tmp
      #  #Do the adjustment
      #  out@x[toAdjust] = out@x[toAdjust]-soupFrac[out@i[toAdjust]+1]*toDo[out@j[toAdjust]+1]
      #  #Only keep updating those that need it
      #  toAdjust = toAdjust[out@x[toAdjust]>0]
      #  out@x[out@x<0]=0
      #  if(verbose>1)
      #    print(quantile(colSums(out)-tgts))
      #  if(max(colSums(out)-tgts)<tol)
      #    break
      #}
      ##This is the clearer but slower version of above
      #out@x = pmax(0,out@x - soupFrac[out@i+1]*expSoupCnts[out@j+1])
      #while(max(colSums(out)-tgts)>tol){
      #  toDo = colSums(out)-tgts
      #  out@x = pmax(0,out@x - soupFrac[out@i+1]*toDo[out@j+1])
      #  print(quantile(colSums(out)-tgts))
      #}
      #Fix the object internals
      w = which(out@x>0)
      out = sparseMatrix(i=out@i[w]+1,
                         j=out@j[w]+1,
                         x=out@x[w],
                         dims=dim(out),
                         dimnames=dimnames(out),
                         giveCsparse=FALSE)
    }else{
      stop("Impossible!")
    }
  }
  #Do stochastic rounding to integers if needed
  if(roundToInt){
    if(verbose>1)
      message("Rounding to integers.")
    #Round to integer below, then probabilistically bump back up
    out@x = floor(out@x)+rbinom(length(out@x),1,out@x-floor(out@x))
  }
  return(out)
}

#' Automatically calculate the contamination fraction
#'
#' The idea of this method is that genes that are highly expressed in the soup and are marker genes for some population can be used to estimate the background contamination.  Marker genes are identified using the tfidf method (see \code{\link{quickMarkers}}).  The contamination fraction is then calculated at the cluster level for each of these genes and clusters are then aggressively pruned to remove those that give implausible estimates.
#'
#' This set of marker genes is filtered to include only those with tf-idf value greater than \code{tfidfMin}.  A higher tf-idf value implies a more specific marker.  Specifically a cut-off t implies that a marker gene has the property that geneFreqGlobal < exp(-t/geneFreqInClust).  See \code{\link{quickMarkers}}.  It may be necessary to decrease this value for data sets with few good markers.
#'
#' This set of marker genes is filtered down to include only the genes that are highly expressed in the soup, controlled by the \code{soupQuantile} parameter.  Genes highly expressed in the soup provide a more precise estimate of the contamination fraction.
#'
#' The pruning of implausible clusters is based on a call to \code{\link{estimateNonExpressingCells}}.  The parameters \code{maximumContamination=max(contaminationRange)} and \code{rhoMaxFDR} are passed to this function.  The defaults set here are calibrated to aggressively prune anything that has even the weakest of evidence that it is genuinely expressed.
#'
#' For each cluster/gene pair the posterior distribution of the contamination fraction is calculated (based on gamma prior, controlled by \code{priorRho} and \code{priorRhoStdDev}).  These posterior distributions are aggregated to produce a final estimate of the contamination fraction. The logic behind this is that estimates from clusters that truly estimate the contamination fraction will cluster around the true value, while erroneous estimates will be spread out across the range (0,1) without a 'preferred value'.  The most probable value of the contamination fraction is then taken as the final global contamination fraction.
#'
#' @param sc The SoupChannel object.
#' @param topMarkers A data.frame giving marker genes.  Must be sorted by decreasing specificity of marker and include a column 'gene' that contains the gene name.  If set to NULL, markers are estimated using \code{\link{quickMarkers}}.
#' @param tfidfMin Minimum value of tfidf to accept for a marker gene.
#' @param soupQuantile Only use genes that are at or above this expression quantile in the soup.  This prevents inaccurate estimates due to using genes with poorly constrained contribution to the background.
#' @param maxMarkers If we have heaps of good markers, keep only the best \code{maxMarkers} of them.
#' @param contaminationRange Vector of length 2 that constrains the contamination fraction to lie within this range.  Must be between 0 and 1.  The high end of this range is passed to \code{\link{estimateNonExpressingCells}} as \code{maximumContamination}.
#' @param rhoMaxFDR False discovery rate passed to \code{\link{estimateNonExpressingCells}}, to test if rho is less than \code{maximumContamination}.
#' @param priorRho Mode of gamma distribution prior on contamination fraction.
#' @param priorRhoStdDev Standard deviation of gamma distribution prior on contamination fraction.
#' @param doPlot Create a plot showing the density of estimates?
#' @param forceAccept Passed to \code{\link{setContaminationFraction}}.  Should we allow very high contamination fractions to be used.
#' @param verbose Be verbose?
#' @seealso quickMarkers
#' @return A modified SoupChannel object where the global contamination rate has been set.  Information about the estimation is also stored in the slot \code{fit}
#'
#' @importFrom stats dgamma qgamma
#' @importFrom graphics abline lines legend plot
autoEstCont = function(sc,topMarkers=NULL,tfidfMin=1.0,soupQuantile=0.90,maxMarkers=100,contaminationRange=c(0.01,0.8),rhoMaxFDR=0.2,priorRho=0.05,priorRhoStdDev=0.10,doPlot=TRUE,forceAccept=FALSE,verbose=TRUE){
    if(!'clusters' %in% colnames(sc$metaData))
        stop("Clustering information must be supplied, run setClusters first.")
    #First collapse by cluster
    s = split(rownames(sc$metaData),sc$metaData$clusters)
    tmp = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
    ssc = sc
    ssc$toc = tmp
    ssc$metaData = data.frame(nUMIs = colSums(tmp),row.names=colnames(tmp))
    ###################
    # Get best markers
    #Get the top N soup Genes
    soupProf = ssc$soupProfile[order(ssc$soupProfile$est,decreasing=TRUE),]
    soupMin = quantile(soupProf$est,soupQuantile)
    #Find or load markers.
    if(is.null(topMarkers)){
        #Refine this to the best markers we can manage
        mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
        #And only the most specific entry for each gene
        mrks = mrks[order(mrks$gene,-mrks$tfidf),]
        mrks = mrks[!duplicated(mrks$gene),]
        #Order by tfidif maxness
        mrks = mrks[order(-mrks$tfidf),]
        #Apply tf-idf cut-off
        mrks = mrks[mrks$tfidf > tfidfMin,]
    }else{
        mrks = topMarkers
    }
    #Filter to include only those that exist in soup
    tgts = rownames(soupProf)[soupProf$est>soupMin]
    #And get the ones that pass our tfidf cut-off
    filtPass = mrks[mrks$gene %in% tgts,]
    tgts = head(filtPass$gene,n=maxMarkers)
    if(verbose)
        message(sprintf("%d genes passed tf-idf cut-off and %d soup quantile filter.  Taking the top %d.",nrow(mrks),nrow(filtPass),length(tgts)))
    #mrks = mrks[mrks$gene %in% tgts,]
    #tgts = head(mrks$gene,nMarks)
    if(length(tgts)==0){
        stop("No plausible marker genes found.  Reduce tfidfMin or soupQuantile")
    }
    if(length(tgts)<10){
        warning("Fewer than 10 marker genes found.  Consider reducing tfidfMin or soupQuantile")
    }
    ############################
    # Get estimates in clusters
    #Get which ones we'd use and where with canonical method
    tmp = as.list(tgts)
    names(tmp) = tgts
    ute = estimateNonExpressingCells(sc,tmp,maximumContamination=max(contaminationRange),FDR=rhoMaxFDR)
    m = rownames(sc$metaData)[match(rownames(ssc$metaData),sc$metaData$clusters)]
    ute = t(ute[m,,drop=FALSE])
    colnames(ute) = rownames(ssc$metaData)
    #Now calculate the observed and expected counts for each cluster for
    expCnts = outer(ssc$soupProfile$est,ssc$metaData$nUMIs)
    rownames(expCnts) = rownames(ssc$soupProfile)
    colnames(expCnts) = rownames(ssc$metaData)
    expCnts = expCnts[tgts,,drop=FALSE]
    #And the observed ones
    obsCnts = ssc$toc[tgts,,drop=FALSE]
    #We're done, but record some extra data for fun and profit
    #Filter out the shite
    #Get the p-value for this being less than 1
    pp = ppois(obsCnts,expCnts*max(contaminationRange),lower.tail=TRUE)
    qq = p.adjust(pp,method='BH')
    qq = matrix(qq,nrow=nrow(pp),ncol=ncol(pp),dimnames=dimnames(pp))
    #Get the cluster level ratio
    rhos = obsCnts/expCnts
    #Index in range
    rhoIdx = t(apply(rhos,1,function(e) order(order(e))))
    #Make a data.frame with everything
    dd = data.frame(gene = rep(rownames(ute),ncol(ute)),
                    passNonExp = as.vector(ute),
                    rhoEst = as.vector(rhos),
                    rhoIdx = as.vector(rhoIdx),
                    obsCnt = as.vector(obsCnts),
                    expCnt = as.vector(expCnts),
                    isExpressedFDR = as.vector(qq)
    )
    dd$geneIdx = match(dd$gene,mrks$gene)
    dd$tfidf = mrks$tfidf[dd$geneIdx]
    dd$soupIdx = match(dd$gene,rownames(soupProf))
    dd$soupExp = soupProf$est[dd$soupIdx]
    dd$useEst = #dd$obsCnt >= minCnts &
        #dd$isExpressedFDR < rhoMaxFDR &
        #dd$rhoIdx <= min(clustPerGene,floor(ncol(rhoIdx)*maxClustFrac)) &
        dd$passNonExp
    #The logic of piling up desity around the true value gets wonky if the number of estimates is low
    if(sum(dd$useEst)<10)
        warning("Fewer than 10 independent estimates, rho estimation is likely to be unstable.  Consider reducing tfidfMin or increasing SoupMin.")
    if(verbose)
        message(sprintf("Using %d independent estimates of rho.",sum(dd$useEst)))
    #Now aggregate the posterior probabilities for the ones we're including
    p.L = function(x,alpha){if(x==0){0}else{qgamma(alpha,x)}}
    p.U = function(x,alpha){qgamma(1-alpha,x+1)}
    alpha=0.95
    alpha=(1-alpha)/2
    dd$rhoHigh=sapply(seq(nrow(dd)),function(e) p.U(dd$obsCnt[e],alpha)/dd$expCnt[e])
    dd$rhoLow=sapply(seq(nrow(dd)),function(e) p.L(dd$obsCnt[e],alpha)/dd$expCnt[e])
    rhoProbes=seq(0,1,.001)
    #Using 95% confidence intervals
    #tmp = sapply(rhoProbes,function(e) {w=which(dd$useEst & dd$tfidf<1.5);sum(e>=dd$rhoLow[w] & e<=dd$rhoHigh[w])/length(w)})
    #Do a posterior estimation instead.  Use gamma prior defined by mode (priorRho) and standard deviation (priorRhoStdDev), which yields a posterior distribution for gamma of the form dgamma(rho,obsCnt+k,scale=theta/(1+theta*expCnts)). Where k and theta are the parameters for prior distribution derived using the above constraints.
    v2 = (priorRhoStdDev/priorRho)**2
    k = 1 +v2**-2/2*(1+sqrt(1+4*v2))
    theta = priorRho/(k-1)
    tmp = sapply(rhoProbes,function(e) {
        tmp = dd[dd$useEst,]
        mean(dgamma(e,k+tmp$obsCnt,scale=theta/(1+theta*tmp$expCnt)))
    })
    #Calculate prior curve
    xx=dgamma(rhoProbes,k,scale=theta)
    #Get estimates
    w = which(rhoProbes>=contaminationRange[1] & rhoProbes<=contaminationRange[2])
    rhoEst = (rhoProbes[w])[which.max(tmp[w])]
    rhoFWHM = range((rhoProbes[w])[which(tmp[w]>=(max(tmp[w])/2))])
    contEst = rhoEst
    if(verbose)
        message(sprintf("Estimated global rho of %.2f",rhoEst))
    ##I think the best way to do this is based on the density.
    #tmp = density(dd$rhoEst[dd$useEst],...)
    #contEst = tmp$x[which.max(tmp$y)]
    if(doPlot){
        plot(rhoProbes,tmp,'l',
             xlim=c(0,1),
             ylim=c(0,max(c(xx,tmp))),
             frame.plot=FALSE,
             xlab='Contamination Fraction',
             ylab='Probability Density')
        #Add prior
        lines(rhoProbes,xx,lty=2)
        abline(v=rhoProbes[which.max(tmp)],col='red')
        legend(x='topright',
               legend=c(sprintf('prior rho %g(+/-%g)',priorRho,priorRhoStdDev),
                        sprintf('post rho %g(%g,%g)',rhoEst,rhoFWHM[1],rhoFWHM[2]),
                        'rho max'),
               lty=c(2,1,1),
               col=c('black','black','red'),
               bty='n')
        #plot(0,
        #     xlim=c(0,1),
        #     ylim=c(0,max(tmp$y)),
        #     type='n',
        #     frame.plot=FALSE,
        #     xlab='Contamination Fraction',
        #     ylab='Density'
        #     )
        #lines(tmp$x,tmp$y)
        #abline(v=contEst,col='red')
    }
    sc$fit = list(dd=dd,
                  priorRho=priorRho,
                  priorRhoStdDev=priorRhoStdDev,
                  posterior = tmp,
                  rhoEst = rhoEst,
                  rhoFWHM = rhoFWHM
    )
    #Set the contamination fraction
    sc = setContaminationFraction(sc,contEst,forceAccept=forceAccept)
    return(sc)
}

#' Calculate the contamination fraction
#'
#' This function computes the contamination fraction using two user-provided bits of information.  Firstly, a list of sets of genes that can be biologically assumed to be absent in at least some cells in your data set.  For example, these might be haemoglobin genes or immunoglobulin genes, which should not be expressed outside of erythroyctes and antibody producing cells respectively.
#'
#' Secondly, this function needs to know which cells definitely do not express the gene sets described above.  Continuing with the haemoglobin example, which are the erythrocytes that are producing haemoglobin mRNAs and which are non-erythrocytes that we can safely assume produce no such genes.  The assumption made is any expression from a gene set in cell marked as a "non-expressor" for that gene set, must be derived from the soup.  Therefore, the level of contamination present can be estimated from the amount of expression of these genes seen in these cells.
#'
#' Most often, the genesets are user supplied based on your knowledge of the experiment and the cells in which they are genuinely expressed is estimated using \code{\link{estimateNonExpressingCells}}.  However, they can also be supplied directly if other information is available.
#'
#' Usually, there is very little variation in the contamination fraction within a channel and very little power to detect the contamination accurately at a single cell level.  As such, the default mode of operation simply estimates one value of the contamination fraction that is applied to all cells in a channel.
#'
#' The global model fits a simple Poisson glm to the aggregated count data across all cells.
#'
#' Finally, note that if you are not able to find a reliable set of genes to use for contamination estimation, or you do not trust the values produced, the contamination fraction can be manually set by the user using \code{\link{setContaminationFraction}}.
#'
#' @param sc A SoupChannel object.
#' @param nonExpressedGeneList A list containing sets of genes which can be assumed to be non-expressed in a subset of cells (see details).
#' @param useToEst A boolean matrix of dimensions ncol(toc) x length(nonExpressedGeneList) indicating which gene-sets should not be assumed to be non-expressed in each cell.  Row names must correspond to the names of \code{nonExpressedGeneList}.  Usually produced by \code{\link{estimateNonExpressingCells}}.
#' @param verbose Print best estimate.
#' @param forceAccept Passed to \code{\link{setContaminationFraction}}.
#' @return A modified version of \code{sc} with estimates of the contamination (rho) added to the metaData table.
#'
#' @importFrom stats coef confint glm poisson quantile
calculateContaminationFraction = function(sc,nonExpressedGeneList,useToEst,verbose=TRUE,forceAccept=FALSE){
    if(!is(sc,'SoupChannel')){
        stop("sc must be a SoupChannel object")
    }
    #Check that you've provided the genes in the right format
    if(!is.list(nonExpressedGeneList))
        stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
    #Check we can estimate
    if(sum(useToEst)==0)
        stop("No cells specified as acceptable for estimation.  useToEst must not be all FALSE")
    #Construct the data.frame to perform inferance on
    df = list()
    for(i in seq_along(nonExpressedGeneList)){
        tgts = nonExpressedGeneList[[i]]
        #Get the soup fraction for this set
        sFrac = sum(sc$soupProfile[tgts,'est'])
        w = rownames(useToEst)[useToEst[,i]]
        if(length(w)>0){
            #Get the counts
            cnts = as.matrix(sc$toc[tgts,w,drop=FALSE])
            df[[i]] = data.frame(row.names=NULL,
                                 cells=colnames(cnts),
                                 geneSet=i,
                                 soupFrac = sFrac,
                                 counts=colSums(cnts),
                                 stringsAsFactors=FALSE)
        }
    }
    df = do.call(rbind,df)
    df$nUMIs = sc$metaData[df$cells,'nUMIs']
    df$expSoupCnts = df$nUMIs * df$soupFrac
    #Make cells a factor, but preserve ordering in sc.  This ensures the the STAN index is the index in sc$metaData
    df$cells = factor(df$cells,levels=rownames(sc$metaData))
    #Fit Poisson GLM with log-link to get global rho
    sc$fit = glm(counts ~ 1,family=poisson(),offset=log(expSoupCnts),data=df)
    #Finally, add it to the meta-data
    sc = setContaminationFraction(sc,exp(coef(sc$fit)),forceAccept=forceAccept)
    tmp=suppressMessages(confint(sc$fit))
    sc$metaData$rhoLow = exp(tmp[1])
    sc$metaData$rhoHigh = exp(tmp[2])
    if(verbose)
        message(sprintf("Estimated global contamination fraction of %0.2f%%",100*exp(coef(sc$fit))))
    return(sc)
}

#' Construct a SoupChannel object
#'
#' Creates a SoupChannel object that contains everything related to the soup estimation of a single channel.
#' @param tod Table of droplets.  A matrix with columns being each droplet and rows each gene.
#' @param toc Table of counts.  Just those columns of \code{tod} that contain cells.
#' @param metaData Meta data pertaining to the cells.  Optional.  Must be a data-frame with rownames equal to column names of \code{toc}.
#' @param calcSoupProfile By default, the soup profile is calculated using \code{\link{estimateSoup}} with default values.  If you want to do something other than the defaults, set this to \code{FALSE} and call \code{\link{estimateSoup}} manually.
#' @param ... Any other named parameters to store.
#' @return A SoupChannel object.
#' @importFrom Matrix colSums
#'
#' @seealso SoupChannelList estimateSoup setSoupProfile setClusters
SoupChannel = function(tod,toc,metaData=NULL,calcSoupProfile=TRUE,...){
    if(!is.null(metaData) & !all(sort(colnames(toc))==sort(rownames(metaData))))
        stop("Rownames of metaData must match column names of table of counts.")
    #Munge everything into a list
    out = list(tod=tod,toc=toc)
    out = c(out,list(...))
    #Create the metadata object
    out$metaData = data.frame(row.names=colnames(toc),
                              nUMIs = colSums(toc)
    )
    #Merge in supplied metaData if it's present
    if(!is.null(metaData)){
        #Drop nUMIs if it exists
        metaData = metaData[,colnames(metaData)!='nUMIs',drop=FALSE]
        out$metaData = cbind(out$metaData,metaData)
    }
    #Get the droplet UMIs as well, as that's a useful thing to have
    out$nDropUMIs = colSums(tod)
    class(out) = c('list','SoupChannel')
    #Estimate the soup
    if(calcSoupProfile)
        out = estimateSoup(out)
    return(out)
}

#' Print method for SoupChannel
#'
#' Prints a summary of a SoupChannel object.
#'
#' @param x A SoupChannel object.
#' @param ... Currently unused.
#' @return Nothing.  Prints message to console.
print.SoupChannel = function(x,...) {
    message(sprintf("Channel with %d genes and %d cells",nrow(x$toc),ncol(x$toc)))
}

#' PBMC 4K meta data
#'
#' Collection of bits of meta data relating to the 10X PBMC 4K data.
#'
#' This data set pertains to the 10X demonstration PBMC 4K data and includes metadata about it in the \code{data.frame} named \code{PBMC_metaData}.
#'
#' \code{PBMC_metaData} was created using Seurat (v2) to calculate a tSNE representation of the data and cluster cells with these commands.
#' \itemize{
#'   \item \code{set.seed(1)}
#'   \item \code{srat = CreateSeuratObject(sc$toc)}
#'   \item \code{srat = NormalizeData(srat)}
#'   \item \code{srat = ScaleData(srat)}
#'   \item \code{srat = FindVariableGenes(srat)}
#'   \item \code{srat = RunPCA(srat,pcs.compute=30)}
#'   \item \code{srat = RunTSNE(srat,dims.use=seq(30))}
#'   \item \code{srat = FindClusters(srat,dims.use=seq(30),resolution=1)}
#'   \item \code{PBMC_metaData = as.data.frame(srat@dr$tsne@cell.embeddings)}
#'   \item \code{colnames(PBMC_metaData) = c('RD1','RD2')}
#'   \item \code{PBMC_metaData$Cluster = factor(srat@meta.data[rownames(PBMC_metaData),'res.1'])}
#'   \item \code{PBMC_metaData$Annotation = factor(c('7'='B','4'='B','1'='T_CD4','2'='T_CD4','3'='T_CD8','5'='T_CD8','6'='NK','8'='NK','0'='MNP','9'='MNP','10'='MNP','11'='?')[as.character(PBMC_metaData$Cluster)])}
#' }
#'
#' @format \code{PBMC_metaData} is a data.frame with 4 columns: RD1, RD2, Cluster, and Annotation.
#' @usage data(PBMC_metaData)
#' @name PBMC_metaData
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
"PBMC_metaData"

#' SoupChannel from PBMC data
#'
#' \code{\link{SoupChannel}} created from 10X demonstration PBMC 4k data.  The cells have been sub-sampled by a factor of 2 to reduce file size of package.
#'
#' \code{PBMC_sc} was created by running the following commands.
#' \itemize{
#'   \item \code{set.seed(1137)}
#'   \item \code{tmpDir = tempdir(check=TRUE)}
#'   \item \code{download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))}
#'   \item \code{download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))}
#'   \item \code{untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)}
#'   \item \code{untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)}
#'   \item \code{library(SoupX)}
#'   \item \code{PBMC_sc = load10X(tmpDir,calcSoupProfile=FALSE)}
#'   \item \code{PBMC_sc = SoupChannel(PBMC_sc$tod,PBMC_sc$toc[,sample(ncol(PBMC_sc$toc),round(ncol(PBMC_sc$toc)*0.5))])}
#' }
#'
#' @format \code{PBMC_sc} is a \code{SoupChannel} object with 33,694 genes and 2,170 cells.
#' @usage data(PBMC_sc)
#' @name PBMC_sc
#' @docType data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
"PBMC_sc"

#' Toy SoupChanel object
#'
#' A \code{\link{SoupChannel}} object created from the toy data used in examples.
#'
#' The toy data is created from a modified version of the extremely reduced \code{Seurat} \code{pbmc_small} dataset.  It includes clusters, tSNE coordinates and a flat estimate of 0.1 contamination.  It includes data for only 226 genes and 62 cells and should not be used for anything other than testing functions as it is not representative of real data in any way.
#'
#' @format \code{scToy} is a \code{SoupChannel} object.
#' @usage data(scToy)
#' @name scToy
#' @docType data
"scToy"

#' Calculate which cells genuinely do not express a particular gene or set of genes
#'
#' Given a list of correlated genes (e.g. Haemoglobin genes, Immunoglobulin genes, etc.), make an attempt to estimate which cells genuinely do not express each of these gene sets in turn.  The central idea is that in cells that are not genuinely producing a class of mRNAs (such as haemoglobin genes), any observed expression of these genes must be due to ambient RNA contamination.  As such, if we can identify these cells, we can use the observed level of expression of these genes to estimate the level of contamination.
#'
#' The ideal way to do this would be to have a prior annotation of your data indicating which cells are (for instance) red blood cells and genuinely expression haemoglobin genes, and which do not and so only express haemoglobin genes due to contamination.  If this is your circumstance, there is no need to run this function, you can instead pass a matrix encoding which cells are haemoglobin expressing and which are not to \code{\link{calculateContaminationFraction}} via the \code{useToEst} parameter.
#'
#' This function will use a conservative approach to excluding cells that it thinks may express one of your gene sets.  This is because falsely including a cell in the set of non-expressing cells may erroneously inflate your estimated contamination, whereas failing to include a genuine non-expressing cell in this set has no significant effect.
#'
#' To this end, this function will exclude any cluster of cells in which any cell is deemed to have genuine expression of a gene set.  Clustering of data is beyond the scope of this package, but can be performed by the user.  In the case of 10X data mapped using cellranger and loaded using \code{\link{load10X}}, the cellranger graph based clustering is automatically loaded and used.
#'
#' To decide if a cell is genuinely expressing a set of genes, a Poisson test is used.  This tests whether the observed expression is greater than \code{maximumContamination} times the expected number of counts for a set of genes, if the cell were assumed to be derived wholly from the background.  This process can be made less conservative (i.e., excluding fewer cells/clusters) by either decreasing the value of the maximum contamination the user believes is plausible (\code{maximumContamination}) or making the significance threshold for the test more strict (by reducing \code{FDR}).
#'
#' @param sc A SoupChannel object.
#' @param nonExpressedGeneList A list containing sets of genes which will be used to estimate the contamination fraction.
#' @param clusters A named vector indicating how to cluster cells.  Names should be cell IDs, values cluster IDs.  If NULL, we will attempt to load it from sc$metaData$clusters.  If set to FALSE, each cell will be considered individually.
#' @param maximumContamination The maximum contamination fraction that you would reasonably expect.  The lower this value is set, the more aggressively cells are excluded from use in estimation.
#' @param FDR A Poisson test is used to identify cells to exclude, this is the false discovery rate it uses.  Higher FDR = more aggressive exclusion.
#' @seealso calculateContaminationFraction plotMarkerMap
#' @return A matrix indicating which cells to be used to estimate contamination for each set of genes.  Typically passed to the \code{useToEst} parameter of \code{\link{calculateContaminationFraction}} or \code{\link{plotMarkerMap}}.
estimateNonExpressingCells = function(sc,nonExpressedGeneList,clusters=NULL,maximumContamination=1.0,FDR=0.05){
    if(!is(sc,'SoupChannel'))
        stop("sc is not a valid SoupChannel object")
    #Get clusters if they exist, if they don't, set to individual cells
    if(is.null(clusters)){
        if('clusters' %in% colnames(sc$metaData)){
            clusters = setNames(as.character(sc$metaData$clusters),rownames(sc$metaData))
        }
    }
    #Using each cell as its own cluster
    if(is.null(clusters) || (length(clusters)==1 && clusters==FALSE)){
        message("No clusters found or supplied, using every cell as its own cluster.")
        clusters = setNames(rownames(sc$metaData),rownames(sc$metaData))
    }
    #Check we have coverage of everything
    if(!all(colnames(sc$toc) %in% names(clusters)))
        stop("Invalid cluster specification.  clusters must be a named vector with all column names in the table of counts appearing.")
    #Convert gene list to genuine list if vector
    if(!is.list(nonExpressedGeneList))
        stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
    #Now work out which clusters to use which genes on
    tgtGns = unique(unlist(nonExpressedGeneList))
    dat = sc$toc[tgtGns,,drop=FALSE]
    cnts = do.call(rbind,lapply(nonExpressedGeneList,function(e) colSums(dat[e,,drop=FALSE])))
    #Work out how many counts we'd expect if the cell were maximally contaminated and all expression came from the contamination
    exp = outer(sc$soupProfile[tgtGns,'est'],sc$metaData$nUMIs*maximumContamination)
    colnames(exp) = colnames(cnts)
    rownames(exp) = tgtGns
    exp = do.call(rbind,lapply(nonExpressedGeneList,function(e) colSums(exp[e,,drop=FALSE])))
    #Identify those cells where a gene is definitely expressed
    s = split(names(clusters),clusters)
    clustExp = ppois(cnts-1,exp,lower.tail=FALSE)
    clustExp = t(apply(clustExp,1,p.adjust,method='BH'))
    clustExp = do.call(rbind,lapply(s,function(e) apply(clustExp[,e,drop=FALSE],1,min)))
    clustExp = clustExp>=FDR
    #Expand it out into a full cell matrix
    clustExp = clustExp[match(clusters,rownames(clustExp)),,drop=FALSE]
    rownames(clustExp) = names(clusters)
    #Check that we actually got som
    if(sum(clustExp)==0){
        warning("No non-expressing cells identified.  Consider setting clusters=FALSE, increasing maximumContamination and/or FDR")
    }
    #A small number found
    if(sum(clustExp)>0 && sum(clustExp)<100){
        warning("Fewer than 100 non-expressing cells identified.  The estimation of the contamination fraction may be inaccurate.  Consider setting clusters=FALSE, increasing maximumContamination and/or FDR")
    }
    return(clustExp)
}

#' Get expression profile of soup
#'
#' This is usually called by \code{\link{SoupChannel}}, rather than directly by the user.  Uses the empty droplets in the range provided to calculate the expression profile of the soup under the assumption that these droplets only contain background.
#'
#' @param sc A \code{SoupChannel} object.
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @param keepDroplets Storing the full table of counts for all droplets uses a lot of space and is really only used to estimate the soup profile.  Therefore, it is dropped after the soup profile has been estimated unless this is set to \code{TRUE}.
#' @return A modified version of \code{sc} with an extra \code{soupProfile} entry containing a data.frame with the soup profile and confidence limits for all genes.
#'
estimateSoup = function(sc,soupRange=c(0,100),keepDroplets=FALSE){
    if(!is(sc,'SoupChannel'))
        stop("sc must be a SoupChannel object.")
    #Estimate the soup
    w = which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
    sc$soupProfile = data.frame(row.names=rownames(sc$tod),
                                est = rowSums(sc$tod[,w,drop=FALSE])/sum(sc$tod[,w]),
                                counts = rowSums(sc$tod[,w,drop=FALSE]))
    #Saves a lot of space if we can drop the droplets now we're done with them
    if(!keepDroplets)
        sc$tod=NULL
    return(sc)
}

#' Load a collection of 10X data-sets
#'
#' Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#'
#' @param dataDir Top level cellranger output directory (the directory that contains the \code{raw_gene_bc_matrices} folder).
#' @param cellIDs Barcodes of droplets that contain cells.  If NULL, use the default cellranger set.
#' @param channelName The name of the channel to store.  If NULL set to either \code{names(dataDir)} or \code{dataDir} is no name is set.
#' @param readArgs A list of extra parameters passed to \code{Seurat::Read10X}.
#' @param includeFeatures If multiple feature types are present, keep only the types mentioned here and collapse to a single matrix.
#' @param verbose Be verbose?
#' @param ... Extra parameters passed to \code{SoupChannel} construction function.
#' @return A SoupChannel object containing the count tables for the 10X dataset.
#' @seealso SoupChannel estimateSoup
#'
#' @importFrom Seurat Read10X
#' @importFrom utils read.csv
load10X = function(dataDir,cellIDs=NULL,channelName=NULL,readArgs=list(),includeFeatures=c('Gene Expression'),verbose=TRUE,...){
    #Work out which version we're dealing with
    isV3 = dir.exists(file.path(dataDir,'raw_feature_bc_matrix'))
    tgt = file.path(dataDir,
                    ifelse(isV3,'raw_feature_bc_matrix','raw_gene_bc_matrices'))
    #Add the reference genome for the non-V3 ones
    if(!isV3)
        tgt = file.path(tgt,list.files(tgt))
    if(verbose)
        message(sprintf("Loading raw count data"))
    dat = do.call(Read10X,c(list(data.dir=tgt),readArgs))
    if(verbose)
        message(sprintf("Loading cell-only count data"))
    if(!is.null(cellIDs)){
        #Do the same sripping that Seurat does on IDs
        if(all(grepl('\\-1$',cellIDs)))
            cellIDs = gsub('\\-1$','',cellIDs)
        #Check we have the IDs
        if(!all(cellIDs %in% colnames(dat)))
            stop("Not all supplied cellIDs found in raw data.")
        datCells = dat[,match(cellIDs,colnames(dat))]
    }else{
        #Work out which ones contain cells
        tgt = file.path(dataDir,
                        ifelse(isV3,'filtered_feature_bc_matrix','filtered_gene_bc_matrices'))
        if(!isV3)
            tgt = file.path(tgt,list.files(tgt))
        datCells = do.call(Read10X,c(list(data.dir=tgt),readArgs))
        #If it's a list of multiple types, have to decide what to include and collapse to one matrix.
        if(is.list(dat)){
            dat = do.call(rbind,dat[includeFeatures])
            datCells = do.call(rbind,datCells[includeFeatures])
        }
    }
    if(verbose)
        message(sprintf("Loading extra analysis data where available"))
    #Get the cluster annotation if available
    mDat = NULL
    tgt = file.path(dataDir,'analysis','clustering','graphclust','clusters.csv')
    if(file.exists(tgt)){
        clusters = read.csv(tgt)
        mDat = data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
    }
    #Add fine grained clusters too if present
    tgt = file.path(dataDir,'analysis','clustering','kmeans_10_clusters','clusters.csv')
    if(file.exists(tgt)){
        clusters = read.csv(tgt)
        mDat$clustersFine = clusters$Cluster
    }
    #Get tSNE if available and point to it
    tgt = file.path(dataDir,'analysis','tsne','2_components','projection.csv')
    if(file.exists(tgt)){
        tsne = read.csv(tgt)
        if(is.null(mDat)){
            mDat = data.frame(tSNE1=tsne$TSNE.1,tSNE2=tsne$TSNE.2,row.names=tsne$Barcode)
        }else{
            mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat),tsne$Barcode)]
            mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat),tsne$Barcode)]
        }
        DR = c('tSNE1','tSNE2')
    }else{
        DR=NULL
    }
    #Ensure rownames of metadata match column names of counts
    if(!is.null(mDat) && any(rownames(mDat)!=colnames(datCells))){
        rownames(mDat) = gsub('-1$','',rownames(mDat))
        if(any(rownames(mDat)!=colnames(datCells)))
            stop("Error matching meta-data to cell names.")
    }
    #Get a name for the channel
    if(is.null(channelName))
        channelName = ifelse(is.null(names(dataDir)),dataDir,names(dataDir))
    channel = SoupChannel(tod = dat,
                          toc = datCells,
                          metaData = mDat,
                          channelName = channelName,
                          dataDir = dataDir,
                          dataType='10X',
                          isV3=isV3,
                          DR=DR,
                          ...
    )
    return(channel)
}

#' SoupX: Profile, quantify and remove ambient RNA expression from droplet based RNA-seq
#'
#' This package implements the method described in REF.  First a few notes about nomenclature:
#' soup - Used a shorthand to refer to the ambient RNA which is contained in the input solution to droplet based RNA-seq experiments and ends up being sequenced along with the cell endogenous RNAs that the experiment is aiming to quantify.
#' channel - This refers to a single run input into a droplet based sequencing platform.  For Chromium 10X 3' sequencing there are currently 8 "channels" per run of the instrument.  Because the profile of the soup depends on the input solution, this is the minimal unit on which the soup should be estimated and subtracted.
#'
#' The essential step in performing background correction is deciding which genes are not expressed in a reasonable fraction of cells.  This is because SoupX estimates the contamination fraction by comparing the expression of these non-expressed genes in droplets containing cells to the soup defined from empty droplets.  For solid tissue, the set of Haemoglobin genes usually works well.  The key properties a gene should have are:
#' - it should be easy to identify when it is truly expressed (i.e., when it's expressed, it should be highly expressed)
#' - it should be highly specific to a certain cell type or group of cell types so that when the expression level is low, you can be confident that the expression is coming from the soup and not a very low level of expression from the cell
#'
#' Spike-in RNAs are the best case scenario.  In the case where you do not have spike-ins and haemoglobin genes are not viable estimators, the user should begin by using the \link{plotMarkerDistribution} function to plot those genes with bi-modal distributions that have a pattern of expression across cells that is consistent with high cell-type specificity.  The user should then select a set of genes that can be used for estimation from this list.  One or two high quality genes is usually sufficient to obtain a good estimate for the average contamination level of a channel.
#'
#' @docType package
#' @name SoupX
#' @import ggplot2
#' @importFrom Matrix colSums rowSums t
#' @importFrom stats lowess spline approx
#' @importFrom grDevices rainbow
#' @importFrom stats approx dpois optimise p.adjust pbinom phyper ppois qbeta
#' @importFrom utils data read.delim setTxtProgressBar txtProgressBar head
#' @importFrom methods as is
NULL
utils::globalVariables(c('RD1','RD2','nUMIs','est','lower','upper','isLogged','MarkerGroup','Values','rho','qVals','logRatio','expSoupCnts','soupProfile'))

#' Plot correlation of expression profiles of soup and aggregated cells
#'
#' Calculates an expression profile by aggregating counts across all cells and plots this (on a log10 scale) against the expression profile of the soup.
#'
#' @param sc A SoupChannel object.
#' @return A ggplot2 object containing the plot.
plotSoupCorrelation = function(sc){
    if(!is(sc,'SoupChannel'))
        stop("sc not a valid SoupChannel object.")
    #Calculate the cell profile
    cellProfile = rowSums(sc$toc)
    cellProfile = (cellProfile/sum(cellProfile))
    df = data.frame(cellProfile,soupProfile=sc$soupProfile$est)
    gg = ggplot(df,aes(log10(cellProfile),log10(soupProfile))) +
        geom_point(alpha=1/3) +
        geom_abline(intercept=0,slope=1) +
        ylab('log10(Soup Expression)')+
        xlab('log10(Aggregate cell Expression)')
    return(gg)
}

#' Plots the distribution of the observed to expected expression for marker genes
#'
#' If each cell were made up purely of background reads, the expression fraction would equal that of the soup.  This plot compares this expectation of pure background to the observed expression fraction in each cell, for each of the groups of genes in \code{nonExpressedGeneList}.  For each group of genes, the distribution of this ratio is plotted across all cells.  A value significantly greater than 1 (0 on log scale) can only be obtained if some of the genes in each group are genuinely expressed by the cell.  That is, the assumption that the cell is pure background does not hold for that gene.
#'
#' This plot is a useful diagnostic for the assumption that a list of genes is non-expressed in most cell types.  For non-expressed cells, the ratio should cluster around the contamination fraction, while for expressed cells it should be elevated.  The most useful non-expressed gene sets are those for which the genes are either strongly expressed, or not expressed at all.  Such groups of genes will show up in this plot as a bimodal distribution, with one mode containing the cells that do not express these genes around the contamination fraction for this channel and another around a value at some value equal to or greater than 0 (1 on non-log scale) for the expressed cells.
#'
#' The red line shows the global estimate of the contamination for each group of markers.  This is usually lower than the low mode of the distribution as there will typically be a non-negligible number of cells with 0 observed counts (and hence -infinity log ratio).
#'
#' If \code{nonExpressedGeneList} is missing, this function will try and find genes that are very specific to different clusters, as these are often the most useful in estimating the contamination fraction.   This is meant only as a heuristic, which can hopefully provide some inspiration as to a class of genes to use to estimation the contamination for your experiment.  Please do **NOT** blindly use the top N genes found in this way to estimate the contamination.  That is, do not feed this list of genes into \code{\link{calculateContaminationFraction}} without any manual consideration or filtering as this *will over-estimate your contamination* (often by a large amount).  For this reason, these gene names are not returned by the function.
#'
#' @param sc A SoupChannel object.
#' @param nonExpressedGeneList Which sets of genes to use to estimate soup (see \code{\link{calculateContaminationFraction}}).
#' @param maxCells Randomly plot only this many cells to prevent over-crowding.
#' @param tfidfMin Minimum specificity cut-off used if finding marker genes (see \code{\link{quickMarkers}}).
#' @param ... Passed to \code{\link{estimateNonExpressingCells}}
#' @importFrom stats setNames
#' @return A ggplot2 object containing the plot.
#'
plotMarkerDistribution = function(sc,nonExpressedGeneList,maxCells=150,tfidfMin=1,...){
    if(!is(sc,'SoupChannel'))
        stop("sc not a valid SoupChannel object.")
    #Get nonExpressedGeneList algorithmically if missing...
    if(missing(nonExpressedGeneList)){
        message("No gene lists provided, attempting to find and plot cluster marker genes.")
        #Get marker genes instead.  Obviously this requires clustering
        if(!'clusters' %in% colnames(sc$metaData))
            stop("Failed as no clusters found!  Clustering must be set via 'setClusters' to find marker genes.")
        #Get top markers
        mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
        #And only the most specific entry for each gene
        mrks = mrks[order(mrks$gene,-mrks$tfidf),]
        mrks = mrks[!duplicated(mrks$gene),]
        #Order by tfidif maxness
        mrks = mrks[order(-mrks$tfidf),]
        #Apply tf-idf cut-off
        mrks = mrks[mrks$tfidf > tfidfMin,]
        message(sprintf("Found %d marker genes",nrow(mrks)))
        #Of the top markers, order by soup
        mrks = mrks[order(sc$soupProfile[mrks$gene,'est'],decreasing=TRUE),]
        #And keep the top 20
        nonExpressedGeneList = mrks$gene[seq(min(nrow(mrks),20))]
        nonExpressedGeneList = setNames(as.list(nonExpressedGeneList),nonExpressedGeneList)
    }
    #Make non-lists into lists
    if(!is.list(nonExpressedGeneList))
        stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
    #Get the non-expressing matrix
    nullMat = estimateNonExpressingCells(sc,nonExpressedGeneList,...)
    #Calculate the ratio to the soup for each marker group in each cell
    obsProfile = t(t(sc$toc)/sc$metaData$nUMIs)
    #Get the ratio
    tst = lapply(nonExpressedGeneList,function(e) colSums(obsProfile[e,,drop=FALSE])/sum(sc$soupProfile[e,'est']))
    #Unlist the thing
    df = data.frame(MarkerGroup = rep(names(tst),lengths(tst)),
                    Barcode=unlist(lapply(tst,names),use.names=FALSE),
                    Values=unlist(tst,use.names=FALSE))
    #Work out which cells to over-plot
    keep = sample(colnames(sc$toc),min(ncol(sc$toc),maxCells))
    #Calculate p-value for each being over some cut-off
    #qVals = do.call(rbind,lapply(nonExpressedGeneList,function(e) p.adjust(pbinom(colSums(sc$toc[e,,drop=FALSE])-1,sc$metaData$nUMIs,minRho*sum(sc$soupProfile[e,'est']),lower.tail=FALSE),method='BH')))
    #df$qVals = qVals[cbind(match(df[,1],rownames(qVals)),match(df[,2],colnames(qVals)))]
    df$nUMIs = sc$metaData[df$Barcode,'nUMIs']
    #Get the expected number of counts
    expCnts = do.call(rbind,lapply(nonExpressedGeneList,function(e) sc$metaData$nUMIs*sum(sc$soupProfile[e,'est'])))
    colnames(expCnts) = rownames(sc$metaData)
    df$expCnts = expCnts[cbind(match(df$MarkerGroup,rownames(expCnts)),match(df$Barcode,colnames(expCnts)))]
    #Set order of marker group as in input
    df$MarkerGroup = factor(df$MarkerGroup,levels=names(nonExpressedGeneList))
    #Add a line estimating the global rho from each group
    #First calculate global rho using nullMat
    globalRhos=c()
    for(i in seq_along(nonExpressedGeneList)){
        if(sum(nullMat[,i])>0){
            tmp = suppressMessages(calculateContaminationFraction(sc,nonExpressedGeneList[i],nullMat[,i,drop=FALSE],forceAccept=TRUE))
            globalRhos = c(globalRhos,tmp$metaData$rho[1])
        }else{
            globalRhos = c(globalRhos,NA)
        }
    }
    names(globalRhos) = names(nonExpressedGeneList)
    globalRhos = data.frame(MarkerGroup = factor(names(globalRhos),levels=names(nonExpressedGeneList)),
                            rho = log10(globalRhos))
    #tmp = df[df$qVals>0.05,]
    #globRhos = sapply(split(tmp,tmp$MarkerGroup),function(e) sum(sc$toc[nonExpressedGeneList[[unique(e$MarkerGroup)]],e$Barcode])/sum(e$nUMIs)/sum(sc$soupProfile[nonExpressedGeneList[[unique(e$MarkerGroup)]],'est']))
    #globRhos = data.frame(MarkerGroup=factor(names(globRhos),levels=names(nonExpressedGeneList)),
    #                      rho = log10(globRhos))
    #Now turn it into a bunch of violin plots
    gg = ggplot(df,aes(MarkerGroup,log10(Values))) +
        geom_violin() +
        geom_jitter(data=df[df$Barcode %in% keep,],aes(size=log10(expCnts)),height=0,width=0.3,alpha=1/2) +
        geom_line(data=globalRhos,aes(MarkerGroup,rho,group=1),colour='red') +
        geom_point(data=globalRhos,aes(MarkerGroup,rho,group=1),colour='red',shape=2) +
        scale_colour_manual(values=c('TRUE'='red','FALSE'='black')) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(colour='expressed\nby cell')+
        ylab('log10(observed/expected)') +
        xlab('Marker group')
    return(gg)
}

#' Plot ratio of observed to expected counts on reduced dimension map
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how likely a set of genes are to be soup derived on that map.  That is, given a set of genes, this function calculates how many counts would be expected if that droplet were nothing but soup and compares that to the observed count.  This is done via a log2 ratio of the two values.  A Poisson test is performed and points that have a statistically significant enrichment over the background (at 5% FDR) are marked.
#'
#' @param sc SoupChannel object.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.  Try and fetch automatically if missing.
#' @param ratLims Truncate log ratios at these values.
#' @param FDR False Discovery Rate for statistical test of enrichment over background.
#' @param useToEst A vector (usually obtained from \code{\link{estimateNonExpressingCells}}), that will be used to mark cells instead of the usual Poisson test.
#' @return A ggplot2 containing the plot.
#'
plotMarkerMap = function(sc,geneSet,DR,ratLims=c(-2,2),FDR=0.05,useToEst=NULL){
    if(!is(sc,'SoupChannel'))
        stop("sc not a valid SoupChannel object.")
    #Try and get DR if missing
    if(missing(DR))
        DR = sc$metaData[,sc$DR]
    #Make sure DR is sensible
    DR = as.data.frame(DR)
    if(ncol(DR)<2)
        stop("Need at least two reduced dimensions.")
    if(!(all(rownames(DR) %in% colnames(sc$toc))))
        stop("rownames of DR need to match column names of sc$toc")
    #Get the ratio of observed to expected
    obs = colSums(sc$toc[geneSet,,drop=FALSE])
    exp = sc$metaData$nUMIs*sum(sc$soupProfile[geneSet,'est'])
    expRatio = obs/exp
    #Add it to the dimensions
    DR$geneRatio = expRatio[rownames(DR)]
    colnames(DR)[1:2] = c('RD1','RD2')
    #Sanitise the values a little
    tgtScale = c(ratLims[1],0,ratLims[2])
    #Rescale to be between zero and 1
    rescaled = (tgtScale-tgtScale[1])/(max(tgtScale)-tgtScale[1])
    DR$logRatio = log10(DR$geneRatio)
    #Keep -Inf as NA as we're not really interested in those that have zero expression
    DR$logRatio[DR$logRatio < ratLims[1]] = ratLims[1]
    DR$logRatio[DR$logRatio > ratLims[2]] = ratLims[2]
    DR$logRatio[DR$geneRatio==0]=NA
    #Calculate the corrected p-value of
    DR$qVals = p.adjust(ppois(obs-1,exp,lower.tail=FALSE),method='BH')[rownames(DR)]
    colVal = 'qVals<FDR'
    if(!is.null(useToEst)){
        DR$useToEst = useToEst
        colVal='useToEst'
    }
    #Create the plot
    gg = ggplot(DR,aes(RD1,RD2)) +
        #Stick NAs underneath
        geom_point(data=DR[is.na(DR$logRatio),],aes_string(colour=colVal),size=0.25) +
        geom_point(data=DR[!is.na(DR$logRatio),],aes_string(fill='logRatio',colour=colVal),size=2.0,shape=21,stroke=0.5) +
        scale_colour_manual(values=c(`FALSE`='black',`TRUE`='#009933'))+
        xlab('ReducedDim1') +
        ylab('ReducedDim2') +
        scale_fill_gradientn(colours = c('blue','white','red'),
                             values = rescaled,
                             guide='colorbar',
                             limits=ratLims
        )
    gg
}

#' Plot maps comparing corrected/raw expression
#'
#' Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how the expression of a geneSet changes after soup correction.
#'
#' @param sc SoupChannel object.
#' @param cleanedMatrix A cleaned matrix to compare against the raw one.  Usually the output of \code{\link{adjustCounts}}.
#' @param geneSet A vector with the names of the genes to aggregate and plot evidence for.
#' @param DR A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.
#' @param dataType How should data be represented.  Binary sets each cell to expressed or not, counts converts everything to counts, soupFrac plots the fraction of the observed counts that are identified as contamination (i.e., (old-new)/old) for each cell and is the default.
#' @param logData Should we log the thing we plot?
#' @return A ggplot2 containing the plot.
#'
plotChangeMap = function(sc,cleanedMatrix,geneSet,DR,dataType=c('soupFrac','binary','counts'),logData=FALSE){
    dataType = match.arg(dataType)
    if(dataType=='binary')
        logData=FALSE
    #Try and get DR if missing
    if(missing(DR))
        DR = sc$metaData[,sc$DR]
    #Make sure DR is sensible
    DR = as.data.frame(DR)
    if(ncol(DR)<2)
        stop("Need at least two reduced dimensions.")
    if(!(all(rownames(DR) %in% colnames(sc$toc))))
        stop("rownames of DR need to match column names of sc$toc")
    colnames(DR)[1:2] = c('RD1','RD2')
    #Simple one panel version
    if(dataType=='soupFrac'){
        df = DR
        old = colSums(sc$toc[geneSet,rownames(df),drop=FALSE])
        new = colSums(cleanedMatrix[geneSet,rownames(df),drop=FALSE])
        relChange = (old-new)/old
        df$old = old
        df$new = new
        df$relChange=relChange
        nom = 'SoupFrac'
        if(logData){
            df$relChange = log10(df$relChange)
            nom = paste0('log10(',nom,')')
            zLims=c(-2,0)
        }else{
            zLims=c(0,1)
        }
        df = df[order(!is.na(df$relChange)),]
        #Truncate to prevent -Inf -> NA coloured dots
        df$relChange[which(df$relChange<zLims[1])] = zLims[1]
        df$relChange[which(df$relChange>zLims[2])] = zLims[2]
        gg = ggplot(df,aes(RD1,RD2)) +
            geom_point(aes(col=relChange),size=0.5) +
            xlab('ReducedDim1') +
            ylab('ReducedDim2') +
            labs(colour=nom) +
            ggtitle('Change in expression due to soup correction')
    }else{
        dfs=list()
        #Get the raw data
        df = DR
        df$correction = 'Uncorrected'
        if(dataType=='binary'){
            df$data = colSums(sc$toc[geneSet,rownames(df),drop=FALSE])>0
        }else if(dataType=='counts'){
            df$data = colSums(sc$toc[geneSet,rownames(df),drop=FALSE])
        }
        if(logData)
            df$data = log10(df$data)
        dfs[['raw']]=df
        #And the corrected expression
        df = DR
        df$correction = 'Corrected'
        if(dataType=='binary'){
            df$data = colSums(cleanedMatrix[geneSet,rownames(df),drop=FALSE])>0
        }else if(dataType=='counts'){
            df$data = colSums(cleanedMatrix[geneSet,rownames(df),drop=FALSE])
        }
        if(logData)
            df$data = log10(df$data)
        dfs[['correctedExpression']]=df
        dfs = do.call(rbind,dfs)
        #Define order to plot
        lvls = c('Uncorrected','Corrected')
        dfs$correction = factor(dfs$correction,levels=lvls[lvls%in%dfs$correction])
        #Stick the NAs at the bottom
        dfs = dfs[order(!is.na(dfs$data)),]
        zLims=c(NA,NA)
        #Now make the plot
        gg = ggplot(dfs,aes(RD1,RD2)) +
            geom_point(aes(colour=data),size=0.5) +
            xlab('ReducedDim1') +
            ylab('ReducedDim2') +
            labs(colour='geneSet') +
            ggtitle('Comparison of before and after correction') +
            facet_wrap(~correction)
    }
    #Make less ugly colour scheme
    if(dataType!='binary')
        gg = gg + scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=zLims)
    gg
}

#' Gets top N markers for each cluster
#'
#' Uses tf-idf ordering to get the top N markers of each cluster.  For each cluster, either the top N or all genes passing the hypergeometric test with the FDR specified, whichever list is smallest.
#'
#' Term Frequency - Inverse Document Frequency is used in natural language processing to identify terms specific to documents.  This function uses the same idea to order genes within a group by how predictive of that group they are.  The main advantage of this is that it is extremely fast and gives reasonable results.
#'
#' To do this, gene expression is binarised in each cell so each cell is either considered to express or not each gene.  That is, we replace the counts with \code{toc > zeroCut}.  The frequency with which a gene is expressed within the target group is compared to the global frequency to calculate the tf-idf score.  We also calculate a multiple hypothesis corrected p-value based on a hypergeometric test, but this is extremely permissive.
#'
#' @param toc Table of counts.  Must be a sparse matrix.
#' @param clusters Vector of length \code{ncol(toc)} giving cluster membership.
#' @param N Number of marker genes to return per cluster.
#' @param FDR False discover rate to use.
#' @param expressCut Value above which a gene is considered expressed.
#' @return data.frame with top N markers (or all that pass the hypergeometric test) and their statistics for each cluster.
#'
#' \dontrun{
#' #Calculate markers from Seurat (v3) object
#' mrks = quickMarkers(srat@assays$RNA@count,srat@active.ident)
#' }
quickMarkers = function(toc,clusters,N=10,FDR=0.01,expressCut=0.9){
    #Convert to the more manipulable format
    toc = as(toc,'dgTMatrix')
    w = which(toc@x>expressCut)
    #Get the counts in each cluster
    clCnts = table(clusters)
    nObs = split(factor(rownames(toc))[toc@i[w]+1],clusters[toc@j[w]+1])
    nObs = sapply(nObs,table)
    #Calculate the observed and total frequency
    nTot = rowSums(nObs)
    tf = t(t(nObs)/as.integer(clCnts[colnames(nObs)]))
    ntf = t(t(nTot - nObs)/as.integer(ncol(toc)-clCnts[colnames(nObs)]))
    idf = log(ncol(toc)/nTot)
    score = tf*idf
    #Calculate p-values
    qvals = lapply(seq_len(ncol(nObs)),function(e)
        p.adjust(phyper(nObs[,e]-1,nTot,ncol(toc)-nTot,clCnts[colnames(nObs)[e]],lower.tail=FALSE),method='BH'))
    qvals = do.call(cbind,qvals)
    colnames(qvals) = colnames(nObs)
    #Get gene frequency of second best
    sndBest = lapply(seq_len(ncol(tf)),function(e) apply(tf[,-e,drop=FALSE],1,max))
    sndBest = do.call(cbind,sndBest)
    colnames(sndBest) = colnames(tf)
    #And the name
    sndBestName = lapply(seq_len(ncol(tf)),function(e) (colnames(tf)[-e])[apply(tf[,-e,drop=FALSE],1,which.max)])
    sndBestName = do.call(cbind,sndBestName)
    colnames(sndBestName) = colnames(tf)
    rownames(sndBestName) = rownames(tf)
    #Now get the top N for each group
    w = lapply(seq_len(ncol(nObs)),function(e){
        o = order(score[,e],decreasing=TRUE)
        if(sum(qvals[,e]<FDR)>=N){
            o[seq(N)]
        }else{
            o[qvals[o,e]<FDR]
        }
    })
    #Now construct the data.frame with everything
    ww = cbind(unlist(w,use.names=FALSE),rep(seq_len(ncol(nObs)),lengths(w)))
    out = data.frame(gene = rownames(nObs)[ww[,1]],
                     cluster = colnames(nObs)[ww[,2]],
                     geneFrequency = tf[ww],
                     geneFrequencyOutsideCluster = ntf[ww],
                     geneFrequencySecondBest = sndBest[ww],
                     geneFrequencyGlobal = nTot[ww[,1]]/ncol(toc),
                     secondBestClusterName = sndBestName[ww],
                     tfidf = score[ww],
                     idf = idf[ww[,1]],
                     qval = qvals[ww],
                     stringsAsFactors=FALSE)
    return(out)
}

#' Set soup profile
#'
#' Manually sets or updates the soup profile for a SoupChannel object.
#'
#' @param sc A SoupChannel object.
#' @param soupProfile A data.frame with columns \code{est} containing the fraction of soup for each gene, \code{counts} containing the total counts for each gene and with row names corresponding to the row names of \code{sc$toc}.
#' @return An updated SoupChannel object with the soup profile set.
#'
setSoupProfile = function(sc,soupProfile){
    if(! 'est' %in% colnames(soupProfile))
        stop("est column missing from soupProfile")
    if(! 'counts' %in% colnames(soupProfile))
        stop("counts column missing from soupProfile")
    if(!all(rownames(soupProfile) %in% rownames(sc$toc))){
        stop("soupProfile invalid.  Not all genes found.")
    }else{
        sc$soupProfile = soupProfile[rownames(sc$toc),]
    }
    return(sc)
}

#' Sets clustering for SoupChannel
#'
#' Adds or updates clustering information to meta-data table in SoupChannel object.
#'
#' @param sc A SoupChannel object.
#' @param clusters A named vector, where entries are the cluster IDs and names are cellIDs.  If no names are provided, the order is assumed to match the order in \code{sc$metaData}.
#' @return An updated SoupChannel object with clustering information stored.
#'
setClusters = function(sc,clusters){
    if(!all(colnames(sc$toc) %in% names(clusters))){
        if(length(clusters)!=nrow(sc$metaData)){
            stop("Invalid cluster specification.  See help.")
        }else{
            sc$metaData$clusters = clusters
        }
    }else{
        sc$metaData$clusters = clusters[rownames(sc$metaData)]
    }
    return(sc)
}

#' Manually set contamination fraction
#'
#' Manually specify the contamination fraction.
#'
#' @param sc A SoupChannel object.
#' @param contFrac The contamination fraction.  Either a constant, in which case the same value is used for all cells, or a named vector, in which case the value is set for each cell.
#' @param forceAccept A warning or error is usually returned for extremely high contamination fractions.  Setting this to TRUE will turn these into messages and proceed.
#' @return A modified SoupChannel object for which the contamination (rho) has been set.
#'
setContaminationFraction = function(sc,contFrac,forceAccept=FALSE){
    #Never want to let this through no matter what
    if(any(contFrac>1))
        stop("Contamination fraction greater than 1 detected.  This is impossible and likely represents a failure in the estimation procedure used.")
    #Let everything else through with a diagnostic message
    if(forceAccept){
        warning = message
        stop = message
    }
    if(any(contFrac>0.5)){
        stop(sprintf("Extremely high contamination estimated (%.2g).  This likely represents a failure in estimating the contamination fraction.  Set forceAccept=TRUE to proceed with this value.",max(contFrac)))
    }else if(any(contFrac>0.3)){
        warning(sprintf("Estimated contamination is very high (%.2g).",max(contFrac)))
    }
    #Now do the setting
    if(length(contFrac)==1){
        sc$metaData$rho=contFrac
    }else{
        if(!all(names(contFrac) %in% rownames(sc$metaData)))
            stop("contFrac must be either of length 1 or a named vector with names matching the rows of sc$metaData")
        sc$metaData$rho[match(names(contFrac),rownames(sc$metaData))] = contFrac
    }
    return(sc)
}

#' Manually set dimension reduction for a channel
#'
#' Manually specify the dimension reduction
#'
#' @param sc A SoupChannel object.
#' @param DR The dimension reduction coordinates (e.g., tSNE).  This must be a data.frame, with two columns giving the two dimension reduction coordinates.  The data.frame must either have row names matching the row names of sc$metaData, or be ordered in the same order as sc$metaData.
#' @param reductName What to name the reduction (defaults to column names provided).
#' @return A modified SoupChannel object for which the dimension reduction has been set.
#'
setDR = function(sc,DR,reductName=NULL){
    #If more than two columns, keep the first two
    if(ncol(DR)>2){
        warning(sprintf("DR has %d columns where 2 were expected.  Using first two.",ncol(DR)))
        DR = DR[,1:2]
    }
    #Check if the rownames match the metadata
    m = match(rownames(sc$metaData),rownames(DR))
    if(any(is.na(m))){
        #Can't use row-names, so the number better match
        if(nrow(DR)!=nrow(sc$metaData)){
            stop(sprintf("Rownames present in metaData not found in DR and row numbers differ (%d in metaData, %d in DR).  Each cell must have a corresponding entry in DR.",nrow(sc$metaData),nrow(DR)))
        }
        m = seq(nrow(DR))
    }
    #Should we change the names?
    if(!is.null(reductName))
        colnames(DR) = paste0(reductName,'_',1:2)
    #Add the entries in
    sc$metaData = cbind(sc$metaData,DR[m,])
    #And point them in the right place
    sc$DR = colnames(DR)
    return(sc)
}

#' Expands soup counts calculated at the cluster level to the cell level
#'
#' Given a clustering of cells and soup counts calculated for each of those clusters, determines a most likely allocation of soup counts at the cell level.
#'
#' @param clustSoupCnts Matrix of genes (rows) by clusters (columns) where counts are number of soup counts for that gene/cluster combination.
#' @param cellObsCnts Matrix of genes (rows) by cells (columns) giving the observed counts
#' @param clusters Mapping from cells to clusters.
#' @param cellWeights Weighting to give to each cell when distributing counts.  This would usually be set to the number of expected soup counts for each cell.
#' @param verbose Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.
#' @return A matrix of genes (rows) by cells (columns) giving the number of soup counts estimated for each cell.  Non-integer values possible.
expandClusters = function(clustSoupCnts,cellObsCnts,clusters,cellWeights,verbose=1){
    ws = cellWeights
    #Do one cluster at a time
    if(verbose>0)
        message(sprintf("Expanding counts from %d clusters to %d cells.",ncol(clustSoupCnts),ncol(cellObsCnts)))
    #lapply instead of loop is a hold-over from when mclapply was an option
    out = lapply(seq(ncol(clustSoupCnts)),
                 function(j) {
                     if(verbose>1)
                         message(sprintf("Expanding cluster %s",colnames(clustSoupCnts)[j]))
                     #Which cells
                     wCells = which(clusters==colnames(clustSoupCnts)[j])
                     #How should they be weighted
                     ww = ws[wCells]/sum(ws[wCells])
                     #What is the limits
                     lims = cellObsCnts[,wCells,drop=FALSE]
                     #And how many soup
                     nSoup = clustSoupCnts[,j]
                     #Create the output object
                     expCnts = as(lims,'dgTMatrix')
                     #Most cases are easily dealt with.  In rough order of frequency.
                     #1. No soup for gene - set to zero
                     #2. All counts for gene are soup - set to lims
                     #3. Not all counts are soup, but every entry is a 0 or 1 so no iteration needed.
                     #4. Some iteration needed.
                     #Deal with case 1
                     expCnts@x[(expCnts@i+1) %in% which(nSoup==0)]=0
                     #Case 2 is dealt with by construction
                     #Get set of genes for cases 3 and 4
                     wGenes = which(nSoup>0 & nSoup<rowSums(lims))
                     #And deal with them as appropriate.  Save time by looking only at non-zero entries
                     w = which((expCnts@i+1) %in% wGenes)
                     w = split(w,expCnts@i[w]+1)
                     tmp = lapply(w,function(e) alloc(nSoup[expCnts@i[e[1]]+1],expCnts@x[e],ww[expCnts@j[e]+1]))
                     expCnts@x[unlist(w,use.names=FALSE)] = unlist(tmp,use.names=FALSE)
                     return(expCnts)
                 })
    out = do.call(cbind,out)
    out = out[,colnames(cellObsCnts)]
    return(out)
}

#' Create Seurat style progress bar
#'
#' Creates progress bar that won't ruin log files and shows progress towards 100%.
#'
#' @param min Minimum value of parameter.
#' @param max Maximum value of parameter.
#' @param ... Passed to \code{\link{txtProgressBar}}
#' @return A txtProgressBar object to use updating progress.
initProgBar = function(min,max,...){
    message('0%   10   20   30   40   50   60   70   80   90   100%')
    message('|----|----|----|----|----|----|----|----|----|----|')
    pb=txtProgressBar(min=min,max=max,style=1,width=51,char='*',...)
    return(pb)
}

#' Allocate values to "buckets" subject to weights and constraints
#'
#' Allocates \code{tgt} of something to \code{length(bucketLims)} different "buckets" subject to the constraint that each bucket has a maximum value of \code{bucketLims} that cannot be exceeded.  By default counts are distributed equally between buckets, but weights can be provided using \code{ws} to have the redistribution prefer certain buckets over others.
#'
#' @param tgt Value to distribute between buckets.
#' @param bucketLims The maximum value that each bucket can take.  Must be a vector of positive values.
#' @param ws Weights to be used for each bucket.  Default value makes all buckets equally likely.
#' @return A vector of the same length as \code{bucketLims} containing values distributed into buckets.
#'
alloc = function(tgt,bucketLims,ws=rep(1/length(bucketLims),length(bucketLims))){
    #Normalise weights
    ws = ws/sum(ws)
    #Save time in line
    if(all(tgt*ws<=bucketLims))
        return(tgt*ws)
    #Need to order things in the order they'll be removed as the tgt increases
    o = order(bucketLims/ws)
    w = ws[o]
    y = bucketLims[o]
    #The formula for number removed at entry i is
    #k_i = \frac{y_i}{w_i} (1- \sum_j=0^{i-1} w_j) + \sum_j=0^{i-1} y_j
    cw = cumsum(c(0,w[-length(w)]))
    cy = cumsum(c(0,y[-length(y)]))
    k = y/w* (1 - cw) + cy
    #Handle zero-weights appropriately
    k[w==0] = Inf
    #Everything that has k<=tgt will be set to y
    b = (k<=tgt)
    #We then need to work out how many counts to distribute we have left over and distribute them according to re-normalised weights
    resid = tgt-sum(y[b])
    w = w/(1-sum(w[b]))
    out = ifelse(b,y,resid*w)
    #Need to reverse sort
    return(out[order(o)])
}
