#' Harmony single cell integration
#'
#' Run Harmony algorithm with Seurat and SingleCellAnalysis pipelines.
#'
#' @param object Pipeline object. Must have PCA computed.
#' @param group.by.vars Which variable(s) to remove (character vector).
#' @param dims.use Which PCA dimensions to use for Harmony. By default, use all
#' @param theta Diversity clustering penalty parameter. Specify for each
#' variable in group.by.vars. Default theta=2. theta=0 does not encourage any
#'  diversity. Larger values of theta result in more diverse clusters.
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#' in group.by.vars. Default lambda=1. Lambda must be strictly positive.
#' Smaller values result in more aggressive correction.
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#' the distance from a cell to cluster centroids. Larger values of sigma result
#'  in cells assigned to more clusters. Smaller values of sigma make soft
#'  kmeans cluster approach hard clustering.
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple
#'  linear regression.
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster.
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param max.iter.cluster Maximum number of rounds to run clustering at each
#' round of Harmony.
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony
#'  Set to -Inf to never stop early.
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step.
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early.
#' @param plot_convergence Whether to print the convergence plot of the
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging.
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to
#'  suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). Cells
#' that have batch variables values matching reference_values will not be moved
#' @param reduction.save Keyword to save Harmony reduction. Useful if you want
#' to try Harmony with multiple parameters and save them as e.g.
#' 'harmony_theta0', 'harmony_theta1', 'harmony_theta2'
#' @param assay.use (Seurat V3 only) Which assay to Harmonize with (RNA
#' by default).
#' @param ... other parameters
#'
#'
#' @rdname RunHarmony
RunHarmony <- function(object, group.by.vars, ...) {
    UseMethod("RunHarmony")
}



#' @rdname RunHarmony
#' @param reduction Name of dimension reduction to use. Default is PCA.
#' @param project.dim Project dimension reduction loadings. Default TRUE.
#' @return Seurat (version 3) object. Harmony dimensions placed into
#' dimensional reduction object harmony. For downstream Seurat analyses,
#' use reduction='harmony'.
RunHarmony.Seurat <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = 'RNA',
  project.dim = TRUE,
  ...
) {
  if (reduction == "pca") {
    tryCatch(
      embedding <- Seurat::Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) {
          message("Harmony needs PCA. Trying to run PCA now.")
        }
        tryCatch(
          object <- Seurat::RunPCA(
            object,
            assay = assay.use, verbose = verbose
          ),
          error = function(e) {
            stop("Harmony needs PCA. Tried to run PCA and failed.")
          }
        )
      }
    )
  } else {
    available.dimreduc <- names(methods::slot(object = object, name = "reductions"))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  }
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(object, group.by.vars)

  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )

  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))

  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(
      embeddings = harmonyEmbed,
      stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
      assay = assay.use,
      key = reduction.save
    )
  })

  object[[reduction.save]] <- harmonydata
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}



#' @rdname RunHarmony
#' @return SingleCellExperiment object. After running RunHarmony, the corrected
#' cell embeddings can be accessed with reducedDim(object, "Harmony").

RunHarmony.SingleCellExperiment <- function(
    object,
    group.by.vars,
    dims.use = NULL,
    theta = NULL,
    lambda = NULL,
    sigma = 0.1,
    nclust = NULL,
    tau = 0,
    block.size = 0.05,
    max.iter.harmony = 10,
    max.iter.cluster = 20,
    epsilon.cluster = 1e-5,
    epsilon.harmony = 1e-4,
    plot_convergence = FALSE,
    verbose = TRUE,
    reference_values = NULL,
    reduction.save = "HARMONY",
    ...
) {

    ## Get PCA embeddings
    if (!"PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("PCA must be computed before running Harmony.")
    }
    pca_embedding <- SingleCellExperiment::reducedDim(object, "PCA")
    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }

    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }
    dims_avail <- seq_len(ncol(pca_embedding))
    if (!all(dims.use %in% dims_avail)) {
        stop("trying to use more dimensions than computed with PCA. Rereun
            PCA with more dimensions or use fewer PCs")
    }

    metavars_df <- SingleCellExperiment::colData(object)
    if (!all(group.by.vars %in% colnames(metavars_df))) {
        stop('Trying to integrate over variables missing in colData')
    }

    harmonyEmbed <- HarmonyMatrix(
        pca_embedding,
        metavars_df,
        group.by.vars,
        FALSE,
        0,
        theta,
        lambda,
        sigma,
        nclust,
        tau,
        block.size,
        max.iter.harmony,
        max.iter.cluster,
        epsilon.cluster,
        epsilon.harmony,
        plot_convergence,
        FALSE,
        verbose,
        reference_values
    )

    rownames(harmonyEmbed) <- row.names(metavars_df)
    colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))
    SingleCellExperiment::reducedDim(object, reduction.save) <- harmonyEmbed

    return(object)
}

#' List of metadata table and scaled PCs matrix
#'
#' @format:
#'   meta_data: data.table of 9478 rows with defining dataset and cell_type
#'   scaled_pcs: data.table of 9478 rows (cells) and 20 columns (PCs)
#'
#' @source \url{support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines"

#' Same as cell_lines but smaller (300 cells).
#'
#' @source \url{support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines_small"

#' Same as cell_lines_small but as Seurat Version 2 object.
#' Expression matrices filled in with dummy values.
#'
#' @source \url{support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines_small_seurat_v2"

#' Same as cell_lines_small but as Seurat Version 3 object.
#' Expression matrices filled in with dummy values.
#'
#' @source \url{support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines_small_seurat_v3"

#' Same as cell_lines_small but as SingleCellExperiment object.
#' Expression matrices filled in with dummy values.
#'
#' @source \url{support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines_small_sce"


#' Harmony: fast, accurate, and robust single cell integration.
#'
#' Algorithm for single cell integration.
#'
#' @section Usage:
#'
#' \enumerate{
#' \item ?HarmonyMatrix to run Harmony on gene expression or PCA
#' embeddings matrix.
#' \item ?RunHarmony to run Harmony on Seurat or SingleCellExperiment objects.
#' }
#' @section Useful links:
#'
#' \enumerate{
#' \item Report bugs at \url{https://github.com/immunogenomics/harmony/issues}
#' \item Read the manuscript
#' \href{https://www.biorxiv.org/content/10.1101/461954v2}{online}.
#' }
#'
#'
#' @name scCancer
#' @docType package
#' @useDynLib scCancer
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp loadModule
#' @importFrom methods new
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data

NULL

#' Main Harmony interface
#'
#' Use this to run the Harmony algorithm on gene expression or PCA matrix.
#'
#' @param data_mat Matrix of normalized gene expession (default) or PCA
#' embeddings (see do_pca).
#' Cells can be rows or columns.
#' @param meta_data Either (1) Dataframe with variables to integrate or (2)
#' vector with labels.
#' @param vars_use If meta_data is dataframe, this defined which variable(s)
#' to remove (character vector).
#' @param do_pca Whether to perform PCA on input matrix.
#' @param npcs If doing PCA on input matrix, number of PCs to compute.
#' @param theta Diversity clustering penalty parameter. Specify for each
#'  variable in vars_use Default theta=2. theta=0 does not encourage any
#'  diversity. Larger values of theta result in more diverse clusters.
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#'  in vars_use.
#' Default lambda=1. Lambda must be strictly positive. Smaller values result
#' in more aggressive correction.
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#'  the distance from a cell to cluster centroids. Larger values of sigma
#'  result in cells assigned to more clusters. Smaller values of sigma make
#'  soft kmeans cluster approach hard clustering.
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple
#' linear regression.
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster.
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param max.iter.cluster Maximum number of rounds to run clustering at each
#' round of Harmony.
#' @param epsilon.cluster Convergence tolerance for clustering round of
#' Harmony. Set to -Inf to never stop early.
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step.
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early.
#' @param plot_convergence Whether to print the convergence plot of the
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging.
#' @param return_object (Advanced Usage) Whether to return the Harmony object
#' or only the corrected PCA embeddings.
#' @param verbose Whether to print progress messages. TRUE to print,
#' FALSE to suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s).
#' Cells that have batch variables values matching reference_values will not
#' be moved.
#' @param cluster_prior (Advanced Usage) Provides user defined clusters for
#' cluster initialization. If the number of provided clusters C is less than K,
#' Harmony will initialize K-C clusters with kmeans. C cannot exceed K.
#'
#' @return By default, matrix with corrected PCA embeddings. If return_object
#' is TRUE, returns the full Harmony object (R6 reference class type).
#'
#'
#'
HarmonyMatrix <- function(
    data_mat, meta_data, vars_use, do_pca = TRUE,
    npcs = 20, theta = NULL, lambda = NULL, sigma = 0.1,
    nclust = NULL, tau = 0, block.size = 0.05,
    max.iter.harmony = 10, max.iter.cluster = 200,
    epsilon.cluster = 1e-5, epsilon.harmony = 1e-4,
    plot_convergence = FALSE, return_object = FALSE,
    verbose = TRUE, reference_values = NULL, cluster_prior = NULL
) {


    ## TODO: check for
    ##    partially observed batch variables (WARNING)
    ##    batch variables with only 1 level (WARNING)
    ##    if lambda given, check correct length
    ##    if theta given, check correct length
    ##    very small batch size and tau=0: suggest tau>0
    ##    is PCA correct?
    if (!(is(meta_data, 'data.frame') | is(meta_data, 'DataFrame'))) {
        #     if (!c('data.frame', '') %in% class(meta_data)) {
        if (length(meta_data) %in% dim(data_mat)) {
            meta_data <- data.frame(batch_variable = meta_data)
            vars_use <- 'batch_variable'
        } else {
            stop('meta_data must be either a data.frame or a vector with batch
                 values for each cell')
        }
        }

    if (is.null(vars_use) | any(!vars_use %in% colnames(meta_data))) {
        msg <- gettextf('must provide variables names (e.g. vars_use=%s)',
                        sQuote('stim'))
        stop(msg)
    }

    if (do_pca) {
        if (ncol(data_mat) != nrow(meta_data)) {
            data_mat <- Matrix::t(data_mat)
        }

        pca_res <- data_mat %>%
            scaleData() %>%
            irlba::prcomp_irlba(n = npcs, retx = TRUE, center = FALSE,
                                scale. = FALSE)
        data_mat <- pca_res$rotation %*% diag(pca_res$sdev)
    }

    N <- nrow(meta_data)
    cells_as_cols <- TRUE
    if (ncol(data_mat) != N) {
        if (nrow(data_mat) == N) {
            data_mat <- t(data_mat)
            cells_as_cols <- FALSE
        } else {
            stop("number of labels do not correspond to number of
                 samples in data matrix")
        }
        }

    if (is.null(nclust)) {
        nclust <- min(round(N / 30), 100)
    }
    if (is.null(theta)) {
        theta <- rep(2, length(vars_use))
    } else if (length(theta) != length(vars_use)) {
        stop('Please specify theta for each variable')
    }
    if (is.null(lambda)) {
        lambda <- rep(1, length(vars_use))
    } else if (length(lambda) != length(vars_use)) {
        stop('Please specify lambda for each variable')
    }
    if (length(sigma) == 1 & nclust > 1) {
        sigma <- rep(sigma, nclust)
    }

    ## Pre-compute some useful statistics
    phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
        t(onehot(meta_data[[var_use]]))
    }))
    N_b <- rowSums(phi)
    Pr_b <- N_b / N
    B_vec <- Reduce(c, lapply(vars_use, function(var_use) {
        length(unique(meta_data[[var_use]]))
    }))
    theta <- Reduce(c, lapply(seq_len(length(B_vec)), function(b)
        rep(theta[b], B_vec[b])))
    theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))

    lambda <- Reduce(c, lapply(seq_len(length(B_vec)), function(b)
        rep(lambda[b], B_vec[b])))
    lambda_mat <- diag(c(0, lambda))

    ## TODO: check that each ref val matches exactly one covariate
    ## TODO: check that you haven't marked all cells as reference!
    if (!is.null(reference_values)) {
        idx <- which(row.names(phi) %in% reference_values)
        cells_ref <- which(colSums(phi[idx, , drop = FALSE] == 1) >= 1)
        b_keep <- which(!row.names(phi) %in% reference_values)
        phi_moe <- phi[b_keep, , drop = FALSE]
        phi_moe[, cells_ref] <- 0

        phi_moe <- rbind(rep(1, N), phi_moe)
        lambda_mat <- lambda_mat[c(1, b_keep + 1), c(1, b_keep + 1)]
    } else {
        phi_moe <- rbind(rep(1, N), phi)
    }

    ## RUN HARMONY
    harmonyObj <- new(harmony, 0) ## 0 is a dummy variable - will change later
    harmonyObj$setup(
        data_mat, phi, phi_moe,
        Pr_b, sigma, theta, max.iter.cluster,epsilon.cluster,
        epsilon.harmony, nclust, tau, block.size, lambda_mat, verbose
    )
    init_cluster(harmonyObj, cluster_prior)
    harmonize(harmonyObj, max.iter.harmony, verbose)
    if (plot_convergence) graphics::plot(HarmonyConvergencePlot(harmonyObj))

    ## Return either the R6 Harmony object or the corrected PCA matrix
    if (return_object) {
        return(harmonyObj)
    } else {
        res <- as.matrix(harmonyObj$Z_corr)
        row.names(res) <- row.names(data_mat)
        colnames(res) <- colnames(data_mat)
        if (!cells_as_cols)
            res <- t(res)
        return(res)
    }
    }

#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom dplyr %>%
#'
#' @usage lhs \%>\% rhs
#' @return return value of rhs function.
NULL


onehot <- function(x) {
    data.frame(x) %>%
        tibble::rowid_to_column("row_id") %>%
        dplyr::mutate(dummy = 1) %>%
        tidyr::spread(x, .data$dummy, fill = 0) %>%
        dplyr::select(-.data$row_id) %>% as.matrix
}

scaleData <- function(A, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- methods::as(A, "dgCMatrix")

    if (margin != 1) A <- t(A)

    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

moe_correct_ridge <- function(harmonyObj) {
    harmonyObj$moe_correct_ridge_cpp()
}


moe_ridge_get_betas <- function(harmonyObj) {
    harmonyObj$moe_ridge_get_betas_cpp()
}


cluster <- function(harmonyObj) {
    if (harmonyObj$ran_init == FALSE) {
        stop('before clustering, run init_cluster')
    }
    harmonyObj$cluster_cpp()
}

harmonize <- function(harmonyObj, iter_harmony, verbose=TRUE) {
    if (iter_harmony < 1) {
        return(0)
    }

    for (iter in seq_len(iter_harmony)) {
        if (verbose) {
            message(gettextf('Harmony %d/%d', iter, iter_harmony))
        }

        # STEP 1: do clustering
        err_status <- cluster(harmonyObj)
        if (err_status == -1) {
            stop('terminated by user')
        } else if (err_status != 0) {
            stop(gettextf('Harmony exited with non-zero exit status: %d',
                          err_status))
        }

        # STEP 2: regress out covariates
        moe_correct_ridge(harmonyObj)

        # STEP 3: check for convergence
        if (harmonyObj$check_convergence(1)) {
            if (verbose) {
                message(gettextf("Harmony converged after %d iterations",
                                 iter))
            }
            return(0)
        }
    }
}


init_cluster <- function(harmonyObj, cluster_prior=NULL) {
    if (harmonyObj$ran_setup == FALSE) {
        stop('before initializing cluster, run setup')
    }
    if (!is.null(cluster_prior)) {
        if (ncol(cluster_prior) != harmonyObj$N) {
            stop('cluster_prior must be defined by N cells')
        }
        if (nrow(cluster_prior) > harmonyObj$K) {
            stop('cluster_prior cannot contain more than K clusters')
        }
        C <- nrow(cluster_prior)
        harmonyObj$Y <- matrix(0, harmonyObj$d, harmonyObj$K)
        harmonyObj$Y[, seq_len(C)] <- compute_Y(harmonyObj$Z_cos,
                                                cluster_prior)
        harmonyObj$R <- matrix(0, harmonyObj$K, harmonyObj$N)
        harmonyObj$R[seq_len(nrow(cluster_prior)), ] <- cluster_prior


        ## if needed, initialize K-C clusters
        if (C < harmonyObj$K) {
            Ynew <- t(stats::kmeans(t(harmonyObj$Z_cos),
                                    centers = harmonyObj$K - C,
                                    iter.max = 25, nstart = 10)$centers)
            harmonyObj$Y[, seq(1+C, harmonyObj$K)] <- Ynew
        }

        harmonyObj$init_cluster_cpp(C)
    } else {
        harmonyObj$Y <- t(stats::kmeans(t(harmonyObj$Z_cos),
                                        centers = harmonyObj$K,
                                        iter.max = 25, nstart = 10)$centers)
        harmonyObj$init_cluster_cpp(0)
    }

}


HarmonyConvergencePlot <- function(
    harmonyObj, round_start=1, round_end=Inf, do_wrap=FALSE
) {
    ## ignore initial value
    ## break down kmeans objective into rounds
    obj_fxn <- data.frame(
        kmeans_idx = Reduce(c, lapply(harmonyObj$kmeans_rounds,
                                      function(rounds) {
                                          seq_len(rounds)
                                      })),
        harmony_idx = Reduce(c, lapply(
            seq_len(length(harmonyObj$kmeans_rounds)),
            function(i) {rep(i, harmonyObj$kmeans_rounds[i])})
        ),
        val = utils::tail(harmonyObj$objective_kmeans, -1)
    ) %>%
        dplyr::filter(.data$harmony_idx >= round_start) %>%
        dplyr::filter(.data$harmony_idx <= round_end) %>%
        tibble::rowid_to_column("idx")


    plt <- obj_fxn %>% ggplot2::ggplot(ggplot2::aes(.data$idx, .data$val,
                                                    col = .data$harmony_idx)) +
        ggplot2::geom_point(shape = 21) +
        ggplot2::labs(y = "Objective Function", x = "Iteration Number")

    if (do_wrap) {
        plt <- plt + ggplot2::facet_grid(.~.data$harmony_idx, scales = 'free',
                                         space = 'free_x')
    }
    return(plt)
}

