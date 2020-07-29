nnmf <- function(
    A, k = 1L, alpha = rep(0,3), beta = rep(0,3), method = c('scd', 'lee'),
    loss = c('mse', 'mkl'), init = NULL, mask = NULL, W.norm = -1L, check.k = TRUE,
    max.iter = 500L, rel.tol = 1e-4, n.threads = 1L, trace = 100/inner.max.iter,
    verbose = 1L, show.warning = TRUE, inner.max.iter = ifelse('mse' == loss, 50L, 1L),
    inner.rel.tol = 1e-9
) {
    method <- match.arg(method);
    loss <- match.arg(loss);
    check.matrix(A, input.name = 'A');
    if (!is.double(A))
        storage.mode(A) <- 'double';
    n <- nrow(A);
    m <- ncol(A);

    init.mask <- reformat.input(init, mask, n, m, k);
    k <- init.mask$K;

    alpha <- c(as.double(alpha), rep(0., 3))[1:3];
    beta <- c(as.double(beta), rep(0., 3))[1:3];
    method.code <- get.method.code(method, loss);

    min.k <- min(dim(A));
    A.isNA <- is.na(A);
    A.anyNA <- any(A.isNA); # anyNA is depreciated in new version of R
    if (A.anyNA) {
        min.k <- min(min.k, ncol(A) - rowSums(A.isNA), nrow(A) - colSums(A.isNA));
    }
    rm(A.isNA);
    if (check.k && k > min.k && all(c(alpha, beta) == 0))
        stop(paste("k larger than", min.k, "is not recommended, unless properly masked or regularized.
                   Set check.k = FALSE if you want to skip this checking."));

    if (n.threads < 0L) n.threads <- 0L; # let openMP decide
    if (is.logical(verbose)) {
        verbose <- as.integer(verbose);
    }
    if (trace <= 0) {
        trace <- 999999L; # only compute error of the 1st and last iteration
    }

    run.time <- system.time(
        out <- c_nnmf(A, as.integer(k),
                      init.mask$Wi, init.mask$Hi, init.mask$Wm, init.mask$Hm,
                      alpha, beta, as.integer(max.iter), as.double(rel.tol),
                      as.integer(n.threads), as.integer(verbose), as.logical(show.warning),
                      as.integer(inner.max.iter), as.double(inner.rel.tol), as.integer(method.code),
                      as.integer(trace))
    );
    names(out) <- c('W', 'H', 'mse', 'mkl', 'target.loss', 'average.epochs', 'n.iteration');
    out$mse <- as.vector(out$mse);
    out$mkl <- as.vector(out$mkl);
    out$target.loss <- as.vector(out$target.loss);
    out$average.epochs <- as.vector(out$average.epochs);

    # add row/col names back
    colnames(out$W) <- colnames(init.mask$Wi);
    rownames(out$H) <- rownames(init.mask$Hi);
    if (!is.null(rownames(A))) rownames(out$W) <- rownames(A);
    if (!is.null(colnames(A))) colnames(out$H) <- colnames(A);
    rm(init.mask);

    if (W.norm > 0) {
        if (is.finite(W.norm)) {
            W.scale <- sapply(out$W, function(x) sum(x^W.norm)^(1./W.norm));
        } else {
            W.scale <- sapply(out$W, max);
        }
        out$W <- out$W %*% diag(1./W.scale);
        out$H <- diag(W.scale) %*% out$H
    }

    out$run.time <- run.time;
    out$options <- list(
        method = method,
        loss = loss,
        alpha = alpha,
        beta = beta,
        init = init,
        mask = mask,
        n.threads = n.threads,
        trace = trace,
        verbose = verbose,
        max.iter = max.iter,
        rel.tol = rel.tol,
        inner.max.iter = inner.max.iter,
        inner.rel.tol = inner.rel.tol
    );
    out$call <- match.call();
    return(structure(out, class = 'nnmf'));
}

get.method.code <- function(method = c('scd', 'lee'), loss = c('mse', 'mkl')) {
    method <- match.arg(method);
    loss <- match.arg(loss);
    code <- 1L;
    if ('mkl' == loss) code <- code + 2L;
    if ('lee' == method) code <- code + 1L;
    return(code);
}

check.matrix <- function(A, dm = NULL, mode = 'numeric', check.na = FALSE, input.name = '', check.negative = FALSE) {
    if (is.null(A)) return(invisible(NULL));
    if (!is.null(dm) && any(dim(A) != dm, na.rm = TRUE))
        stop(sprintf("Dimension of matrix %s is expected to be (%d, %d), but got (%d, %d)", input.name, nrow(A), ncol(A), dm[1], dm[2]));
    if (mode(A) != mode) stop(sprintf("Matrix %s must be %s.", input.name, mode));
    if (check.negative && any(A[!is.na(A)] < 0)) stop(sprintf("Matrix %s must be non-negative.", input.name));
    if (check.na && any(is.na(A))) stop(sprintf("Matrix %s contains missing values.", input.name));
}

reformat.input <- function(init, mask, n, m, k) {
    if (is.null(mask)) mask <- list();
    if (is.null(init)) init <- list();
    stopifnot(is.list(mask));
    stopifnot(is.list(init));

    known.W <- !is.null(init[['W0']]);
    known.H <- !is.null(init[['H0']]);
    kW0 <- kH0 <- 0;

    is.empty <- function(x) 0 == length(x);

    if (known.W) {
        if(!is.matrix(init[['W0']]))
            init[['W0']] <- as.matrix(init[['W0']]);
        kW0 <- ncol(init[['W0']]);
        mask[['W0']] <- matrix(TRUE, n, kW0);
    }
    else {
        mask[['W0']] <- NULL;
        mask[['H1']] <- NULL;
        init[['H1']] <- NULL;
    }

    if (known.H) {
        if(!is.matrix(init$H0))
            init[['H0']] <- as.matrix(init[['H0']]);
        kH0 <- nrow(init[['H0']]);
        mask[['H0']] <- matrix(TRUE, kH0, m);
    }
    else {
        mask[['H0']] <- NULL;
        mask[['W1']] <- NULL;
        init[['W1']] <- NULL;
    }

    K <- k + kW0 + kH0;

    ew <- !all(sapply(mask[c('W', 'W0', 'W1')], is.empty));
    eh <- !all(sapply(mask[c('H', 'H0', 'H1')], is.empty));
    dim.mask <- list(
        'W' = c(n, k*ew), 'W0' = c(n, kW0*ew), 'W1' = c(n, kH0*ew),
        'H' = c(k*eh, m), 'H1' = c(kW0*eh, m), 'H0' = c(kH0*eh, m)
    );

    for (mat in c('W', 'W0', 'W1', 'H' ,'H0', 'H1')) {
        check.matrix(mask[[mat]], dim.mask[[mat]], 'logical', TRUE, paste0('mask$', mat));
        if (is.empty(mask[[mat]]))
            mask[[mat]] <- matrix(FALSE, dim.mask[[mat]][[1]], dim.mask[[mat]][[2]]);
    }

    ew <- !all(sapply(init[c('W', 'W0', 'W1')], is.empty));
    eh <- !all(sapply(init[c('H', 'H0', 'H1')], is.empty));
    dim.init <- list(
        'W' = c(n, k*ew), 'W0' = c(n, kW0*ew), 'W1' = c(n, kH0*ew),
        'H' = c(k*eh, m), 'H1' = c(kW0*eh, m), 'H0' = c(kH0*eh, m)
    );
    for (mat in c('W', 'W0', 'W1', 'H' ,'H0', 'H1')) {
        check.matrix(init[[mat]], dim.init[[mat]], 'numeric', TRUE, paste0('init$', mat));
        if (is.empty(init[[mat]])) {
            init[[mat]] <- matrix(
                runif(prod(dim.init[[mat]])),
                dim.init[[mat]][[1]],
                dim.init[[mat]][[2]]
            );
            #init[[mat]][mask[[mat]]] <- 0;
        }
        if (!is.double(init[[mat]]))
            storage.mode(init[[mat]]) <- 'double';
    }

    return(
        list(
            Wm = do.call(cbind, mask[c('W', 'W0', 'W1')]),
            Hm = do.call(rbind, mask[c('H', 'H1', 'H0')]),
            Wi = do.call(cbind, init[c('W', 'W0', 'W1')]),
            Hi = do.call(rbind, init[c('H', 'H1', 'H0')]),
            kW0 = kW0,
            kH0 = kH0,
            K = K
        ));
}

