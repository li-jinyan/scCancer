#include "harmony.h"
#include "utils.h"
#include "nnlm.h"

//[[Rcpp::export]]
Rcpp::List c_nnmf(const arma::mat & A, const unsigned int k, arma::mat W, arma::mat H, arma::umat Wm, arma::umat Hm,
                  const arma::vec & alpha, const arma::vec & beta, const unsigned int max_iter, const double rel_tol,
                  const int n_threads, const int verbose, const bool show_warning, const unsigned int inner_max_iter,
                  const double inner_rel_tol, const int method, unsigned int trace)
{

    unsigned int n = A.n_rows;
    unsigned int m = A.n_cols;
    //int k = H.n_rows; // decomposition rank k
    unsigned int N_non_missing = n*m;

    if (trace < 1) trace = 1;
    unsigned int err_len = (unsigned int)std::ceil(double(max_iter)/double(trace)) + 1;
    vec mse_err(err_len), mkl_err(err_len), terr(err_len), ave_epoch(err_len);

    // check progression
    bool show_progress = false;
    if (verbose == 1) show_progress = true;
    Progress prgrss(max_iter, show_progress);

    double rel_err = rel_tol + 1;
    double terr_last = 1e99;
    uvec non_missing;
    bool any_missing = !A.is_finite();
    if (any_missing)
    {
        non_missing = find_finite(A);
        N_non_missing = non_missing.n_elem;
        mkl_err.fill(mean((A.elem(non_missing)+TINY_NUM) % log(A.elem(non_missing)+TINY_NUM) - A.elem(non_missing)));
    }
    else
        mkl_err.fill(mean(mean((A+TINY_NUM) % log(A+TINY_NUM) - A))); // fixed part in KL-dist, mean(A log(A) - A)

    if (Wm.empty())
        Wm.resize(0, n);
    else
        inplace_trans(Wm);
    if (Hm.empty())
        Hm.resize(0, m);

    if (W.empty())
    {
        W.randu(k, n);
        W *= 0.01;
        if (!Wm.empty())
            W.elem(find(Wm > 0)).fill(0.0);
    }
    else
        inplace_trans(W);

    if (H.empty())
    {
        H.randu(k, m);
        H *= 0.01;
        if (!Hm.empty())
            H.elem(find(Hm > 0)).fill(0.0);
    }

    if (verbose == 2)
    {
        Rprintf("\n%10s | %10s | %10s | %10s | %10s\n", "Iteration", "MSE", "MKL", "Target", "Rel. Err.");
        Rprintf("--------------------------------------------------------------\n");
    }

    int total_raw_iter = 0;
    unsigned int i = 0;
    unsigned int i_e = 0; // index for error checking
    for(; i < max_iter && std::abs(rel_err) > rel_tol; i++)
    {
        Rcpp::checkUserInterrupt();
        prgrss.increment();

        if (any_missing)
        {
            // update W
            total_raw_iter += update_with_missing(W, H, A.t(), Wm, alpha, inner_max_iter, inner_rel_tol, n_threads, method);
            // update H
            total_raw_iter += update_with_missing(H, W, A, Hm, beta, inner_max_iter, inner_rel_tol, n_threads, method);

            if (i % trace == 0)
            {
                const mat & Ahat = W.t()*H;
                mse_err(i_e) = mean(square((A - Ahat).eval().elem(non_missing)));
                mkl_err(i_e) += mean((-(A+TINY_NUM) % log(Ahat+TINY_NUM) + Ahat).eval().elem(non_missing));
            }
        }
        else
        {
            // update W
            total_raw_iter += update(W, H, A.t(), Wm, alpha, inner_max_iter, inner_rel_tol, n_threads, method);
            // update H
            total_raw_iter += update(H, W, A, Hm, beta, inner_max_iter, inner_rel_tol, n_threads, method);

            if (i % trace == 0)
            {
                const mat & Ahat = W.t()*H;
                mse_err(i_e) = mean(mean(square((A - Ahat))));
                mkl_err(i_e) += mean(mean(-(A+TINY_NUM) % log(Ahat+TINY_NUM) + Ahat));
            }
        }

        if (i % trace == 0)
        {
            ave_epoch(i_e) = double(total_raw_iter)/(n+m);
            if (method < 3) // mse based
                terr(i_e) = 0.5*mse_err(i_e);
            else // KL based
                terr(i_e) = mkl_err(i_e);

            add_penalty(i_e, terr, W, H, N_non_missing, alpha, beta);

            rel_err = 2*(terr_last - terr(i_e)) / (terr_last + terr(i_e) + TINY_NUM );
            terr_last = terr(i_e);
            if (verbose == 2)
                Rprintf("%10d | %10.4f | %10.4f | %10.4f | %10.g\n", i+1, mse_err(i_e), mkl_err(i_e), terr(i_e), rel_err);

            total_raw_iter = 0; // reset to 0
            ++i_e;
        }
    }

    // compute error of the last iteration
    if ((i-1) % trace != 0)
    {
        if (any_missing)
        {
            const mat & Ahat = W.t()*H;
            mse_err(i_e) = mean(square((A - Ahat).eval().elem(non_missing)));
            mkl_err(i_e) += mean((-(A+TINY_NUM) % log(Ahat+TINY_NUM) + Ahat).eval().elem(non_missing));
        }
        else
        {
            const mat & Ahat = W.t()*H;
            mse_err(i_e) = mean(mean(square((A - Ahat))));
            mkl_err(i_e) += mean(mean(-(A+TINY_NUM) % log(Ahat+TINY_NUM) + Ahat));
        }

        ave_epoch(i_e) = double(total_raw_iter)/(n+m);
        if (method < 3) // mse based
            terr(i_e) = 0.5*mse_err(i_e);
        else // KL based
            terr(i_e) = mkl_err(i_e);
        add_penalty(i_e, terr, W, H, N_non_missing, alpha, beta);

        rel_err = 2*(terr_last - terr(i_e)) / (terr_last + terr(i_e) + TINY_NUM );
        terr_last = terr(i_e);
        if (verbose == 2)
            Rprintf("%10d | %10.4f | %10.4f | %10.4f | %10.g\n", i+1, mse_err(i_e), mkl_err(i_e), terr(i_e), rel_err);

        ++i_e;
    }

    if (verbose == 2)
    {
        Rprintf("--------------------------------------------------------------\n");
        Rprintf("%10s | %10s | %10s | %10s | %10s\n\n", "Iteration", "MSE", "MKL", "Target", "Rel. Err.");
    }

    if (i_e < err_len)
    {
        mse_err.resize(i_e);
        mkl_err.resize(i_e);
        terr.resize(i_e);
        ave_epoch.resize(i_e);
    }

    if (show_warning && rel_err > rel_tol)
        Rcpp::warning("Target tolerance not reached. Try a larger max.iter.");

    return Rcpp::List::create(
        Rcpp::Named("W") = W.t(),
        Rcpp::Named("H") = H,
        Rcpp::Named("mse_error") = mse_err,
        Rcpp::Named("mkl_error") = mkl_err,
        Rcpp::Named("target_error") = terr,
        Rcpp::Named("average_epoch") = ave_epoch,
        Rcpp::Named("n_iteration") = i
    );
}

int update(mat & H, const mat & Wt, const mat & A, const umat & mask,
           const vec & beta, unsigned int max_iter, double rel_tol, int n_threads, int method)
{
    // A = W H, solve H
    // No missing in A, Wt = W^T
    // method: 1 = scd, 2 = lee_ls, 3 = scd_kl, 4 = lee_kl

    unsigned int m = A.n_cols;
    int total_raw_iter = 0;

    if (n_threads < 0) n_threads = 0;
    bool is_masked = !mask.empty();
    mat WtW;
    vec mu, sumW;
    if (method == 1 || method == 2)
    {
        WtW = Wt*Wt.t();
        if (beta(0) != beta(1))
            WtW.diag() += beta(0) - beta(1);
        if (beta(1) != 0)
            WtW += beta(1);
        WtW.diag() += TINY_NUM; // for stability: avoid divided by 0 in scd_ls, scd_kl
    }
    else
        sumW = sum(Wt, 1);

#pragma omp parallel for num_threads(n_threads) schedule(dynamic) private(mu)
    for (unsigned int j = 0; j < m; j++) // by columns of H
    {
        // break if all entries of col_j are masked
        if (is_masked && arma::all(mask.col(j)))
            continue;

        int iter = 0;
        if (method == 1)
        {
            mu = WtW*H.col(j) - Wt*A.col(j);
            if (beta(2) != 0)
                mu += beta(2);
            iter = scd_ls_update(H.col(j), WtW, mu, mask.col(j), max_iter, rel_tol);
        }
        else if (method == 2)
            iter = lee_ls_update(H.col(j), WtW, Wt*A.col(j), beta(2), mask.col(j), max_iter, rel_tol);
        else if (method == 3)
            iter = scd_kl_update(H.col(j), Wt, A.col(j), sumW, mask.col(j), beta, max_iter, rel_tol);
        else if (method == 4)
            iter = lee_kl_update(H.col(j), Wt, A.col(j), sumW, mask.col(j), beta, max_iter, rel_tol);

#pragma omp critical
        total_raw_iter += iter;
    }
    return total_raw_iter;
}


int update_with_missing(mat & H, const mat & Wt, const mat & A, const umat & mask,
                        const vec & beta, unsigned int max_iter, double rel_tol, int n_threads, int method)
{
    // A = W H, solve H
    // With missings in A, Wt = W^T
    // method: 1 = scd, 2 = lee_ls, 3 = scd_kl, 4 = lee_kl

    unsigned int n = A.n_rows, m = A.n_cols;
    unsigned int total_raw_iter = 0;

    if (n_threads < 0) n_threads = 0;
    bool is_masked = mask.n_elem > 0;
    mat WtW;
    vec mu;

#pragma omp parallel for num_threads(n_threads) schedule(dynamic) private(WtW, mu)
    for (unsigned int j = 0; j < m; j++) // by columns of H
    {
        // break if all entries of col_j are masked
        if (is_masked && arma::all(mask.col(j)))
            continue;

        bool any_missing = !is_finite(A.col(j));
        uvec non_missing;
        if (any_missing)
            non_missing = find_finite(A.col(j));

        if (method == 1 || method == 2)
        {
            if (any_missing)
            {
                //non_missing.print("Non Missing");
                WtW = Wt.cols(non_missing)*Wt.cols(non_missing).t();
                mu = Wt.cols(non_missing) * A.elem(j*n + non_missing);
            }
            else
            {
                WtW = Wt*Wt.t();
                mu = Wt*A.col(j);
            }
            if (beta(0) != beta(1))
                WtW.diag() += beta(0) - beta(1);
            if (beta(1) != 0)
                WtW += beta(1);
        }

        int iter = 0;
        if (method == 1)
        {
            mu = WtW*H.col(j)-mu;
            if (beta(2) != 0)
                mu += beta(2);
            iter = scd_ls_update(H.col(j), WtW, mu, mask.col(j), max_iter, rel_tol);
        }
        else if (method == 2)
        {
            iter = lee_ls_update(H.col(j), WtW, mu, beta(2), mask.col(j), max_iter, rel_tol);
        }
        else if (method == 3)
        {
            if (any_missing)
                iter = scd_kl_update(H.col(j), Wt.cols(non_missing), A.elem(j*n + non_missing),
                                     sum(Wt.cols(non_missing), 1), mask.col(j), beta, max_iter, rel_tol);
            else
                iter = scd_kl_update(H.col(j), Wt, A.col(j), sum(Wt, 1), mask.col(j), beta, max_iter, rel_tol);
        }
        else if (method == 4)
        {
            if (any_missing)
                iter = lee_kl_update(H.col(j), Wt.cols(non_missing), A.elem(j*n + non_missing),
                                     sum(Wt.cols(non_missing), 1), mask.col(j), beta, max_iter, rel_tol);
            else
                iter = lee_kl_update(H.col(j), Wt, A.col(j), sum(Wt, 1), mask.col(j), beta, max_iter, rel_tol);
        }

#pragma omp critical
        total_raw_iter += iter;
    }
    return total_raw_iter;
}

int scd_ls_update(subview_col<double> Hj, const mat & WtW, vec & mu, const subview_col<uword> mask, const unsigned int & max_iter, const double & rel_tol)
{
    // Problem:  Aj = W * Hj
    // Method: sequential coordinate-wise descent when loss function = square error
    // WtW = W^T W
    // WtAj = W^T Aj
    // beta3: L1 regularization
    // mask: skip updating

    double tmp;
    double etmp = 0;
    double rel_err = 1 + rel_tol;
    bool is_masked = mask.n_elem > 0;

    unsigned int t = 0;
    for (; t < max_iter && rel_err > rel_tol; t++)
    {
        rel_err = 0;
        for (unsigned int k = 0; k < WtW.n_cols; k++)
        {
            if (is_masked && mask(k) > 0) continue;
            tmp = Hj(k) - mu(k) / WtW(k,k);
            if (tmp < 0) tmp = 0;
            if (tmp != Hj(k))
                mu += (tmp - Hj(k)) * WtW.col(k);
            else
                continue;
            etmp = 2*std::abs(Hj(k)-tmp) / (tmp+Hj(k)+TINY_NUM);
            if (etmp > rel_err)
                rel_err = etmp;
            Hj(k) = tmp;
        }
    }
    return int(t);
}


int lee_ls_update(subview_col<double> Hj, const mat & WtW, const vec & WtAj, const double & beta3,
                  const subview_col<uword> mask, const unsigned int & max_iter, const double & rel_tol)
{
    // Problem:  Aj = W * Hj
    // Method: Lee's multiplicative update when loss function = square error
    // WtW = W^T W
    // WtAj = W^T Aj
    // beta3: L1 regularization
    // mask: skip updating

    double tmp;
    double rel_err = rel_tol + 1;
    bool is_masked = mask.n_elem > 0;
    unsigned int t = 0;
    for (; t < max_iter && rel_err > rel_tol; t++)
    {
        rel_err = 0;
        for (unsigned int k = 0; k < WtW.n_cols; k++)
        {
            if (is_masked && mask(k) > 0) continue;
            tmp = dot(WtW.col(k), Hj) + beta3;
            tmp = WtAj(k) / (tmp+TINY_NUM);
            Hj(k) *= tmp;
            tmp = 2*std::abs(tmp-1)/(tmp+1);
            if (tmp > rel_err) rel_err = tmp;
        }
    }
    return int(t);
}


int scd_kl_update(subview_col<double> Hj, const mat & Wt, const vec & Aj, const vec & sumW, const subview_col<uword> mask,
                  const vec & beta, const unsigned int & max_iter, const double & rel_tol)
{
    // Problem:  Aj = W * Hj
    // Method: Sequentially minimize KL distance using quadratic approximation
    // Wt = W^T
    // sumW = column sum of W
    // mask: skip updating
    // beta: a vector of 3, for L2, angle, L1 regularization

    double sumHj = sum(Hj);
    vec Ajt = Wt.t()*Hj;
    vec mu;
    double a; // 2nd-order-derivative
    double b; // 1st-order-derivative
    double tmp, etmp;
    double rel_err = 1 + rel_tol;
    bool is_masked = mask.n_elem > 0;

    unsigned int t = 0;
    for (; t < max_iter && rel_err > rel_tol; t++)
    {
        rel_err = 0;
        for (unsigned int k = 0; k < Wt.n_rows; k++)
        {
            if (is_masked && mask(k) > 0) continue;
            mu = Wt.row(k).t()/(Ajt + TINY_NUM);
            a = dot(Aj, square(mu));
            b = dot(Aj, mu) - sumW(k); // 0.5*ax^2 - bx
            a += beta(0);
            b += a*Hj(k) - beta(2) - beta(1)*(sumHj - Hj(k));
            tmp = b/(a+TINY_NUM);
            if (tmp < 0) tmp = 0;
            if (tmp != Hj(k))
            {
                Ajt += (tmp - Hj(k)) * Wt.row(k).t();
                etmp = 2*std::abs(Hj(k)-tmp) / (tmp+Hj(k) + TINY_NUM);
                if (etmp > rel_err)
                    rel_err = etmp;
                sumHj += tmp - Hj(k);
                Hj(k) = tmp;
            }
        }
    }
    return int(t);
}


int lee_kl_update(subview_col<double> Hj, const mat & Wt, const vec & Aj, const vec & sumW,
                  const subview_col<uword> mask, const vec & beta, const unsigned int & max_iter, const double & rel_tol)
{
    // Problem:  Aj = W * Hj
    // Method: Lee's multiplicative updating when loss = KL divergence
    // Wt = W^T
    // sumW = column sum of W
    // mask: skip updating
    // beta: a vector of 3, for L2, angle, L1 regularization

    double sumHj = sum(Hj);
    double rel_err = rel_tol + 1;
    double tmp;
    bool is_masked = mask.n_elem > 0;
    vec wh = Wt.t()*Hj;
    unsigned int t = 0;
    for (; t < max_iter && rel_err > rel_tol; t++)
    {
        rel_err = 0;
        for (unsigned int k = 0; k < Wt.n_rows; k++)
        {
            if (is_masked && mask(k) > 0) continue;
            tmp = as_scalar(Wt.row(k)*(Aj/(wh+TINY_NUM)));
            tmp /= (sumW(k) + beta(0)*Hj(k) + beta(1)*(sumHj-Hj(k)) + beta(2));
            wh += (tmp-1)*Hj(k) * Wt.row(k).t();
            sumHj += (tmp-1)*Hj(k);
            Hj(k) *= tmp;
            tmp = 2*std::abs(tmp-1)/(tmp+1);
            if (tmp > rel_err) rel_err = tmp;
        }
    }
    return int(t);
}

// add_penalty to the target error 'terr'
void add_penalty(const unsigned int & i_e, vec & terr, const mat & W, const mat & H,
                 const unsigned int & N_non_missing, const vec & alpha, const vec & beta)
{
    // add penalty term back to the loss function (terr)
    if (alpha(0) != alpha(1))
        terr(i_e) += 0.5*(alpha(0)-alpha(1))*accu(square(W))/N_non_missing;
    if (beta(0) != beta(1))
        terr(i_e) += 0.5*(beta(0)-beta(1))*accu(square(H))/N_non_missing;
    if (alpha(1) != 0)
        terr(i_e) += 0.5*alpha(1)*accu(W*W.t())/N_non_missing;
    if (beta(1) != 0)
        terr(i_e) += 0.5*beta(1)*accu(H*H.t())/N_non_missing;
    if (alpha(2) != 0)
        terr(i_e) += alpha(2)*accu(W)/N_non_missing;
    if (beta(2) != 0)
        terr(i_e) += beta(2)*accu(H)/N_non_missing;
}


// NOTE: This is a dummy constructor, needed by Rcpp
harmony::harmony(int __K): K(__K) {}



void harmony::setup(MATTYPE& __Z, MATTYPE& __Phi, MATTYPE& __Phi_moe, VECTYPE __Pr_b,
                    VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans,
                    float __epsilon_kmeans, float __epsilon_harmony,
                    int __K, float tau, float __block_size,
                    MATTYPE __lambda, bool __verbose) {

  Z_corr = MATTYPE(__Z);
  Z_orig = MATTYPE(__Z);
  Z_cos = MATTYPE(Z_orig);
  cosine_normalize(Z_cos, 0, true); // normalize columns

  Phi = __Phi;
  Phi_moe = __Phi_moe;
  N = Z_corr.n_cols;
  Pr_b = __Pr_b;
  B = Phi.n_rows;
  d = Z_corr.n_rows;
  window_size = 3;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_harmony = __epsilon_harmony;


  lambda = __lambda;
  sigma = __sigma;
  sigma_prior = __sigma;
  block_size = __block_size;
  K = __K;
  max_iter_kmeans = __max_iter_kmeans;
  verbose = __verbose;

  theta = __theta;
  allocate_buffers();
  ran_setup = true;
}


void harmony::allocate_buffers() {
  _scale_dist = zeros<MATTYPE>(K, N);
  dist_mat = zeros<MATTYPE>(K, N);
  O = zeros<MATTYPE>(K, B);
  E = zeros<MATTYPE>(K, B);
  W = zeros<MATTYPE>(B + 1, d);
  Phi_Rk = zeros<MATTYPE>(B + 1, N);
}




void harmony::init_cluster_cpp(unsigned C) {
  // kmeans is called outside, in the R function
  cosine_normalize(Y, 0, false); // normalize columns

  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos);

  // if C > 0, only initialize the clusters not set by the user
  // with cluster_prior
  if (C > 0 && C < K) {
      MATTYPE Rtmp = -dist_mat.rows(C, K-1);
      Rtmp.each_col() /= sigma.rows(C, K-1);
      Rtmp.each_row() -= max(Rtmp, 0);
      Rtmp = exp(Rtmp);
      Rtmp.each_row() /= sum(Rtmp, 0);
      R.rows(C, K-1) = Rtmp;
  } else {
      R = -dist_mat;
      R.each_col() /= sigma;
      R.each_row() -= max(R, 0);
      R = exp(R);
      R.each_row() /= sum(R, 0);
  }

  // (3) BATCH DIVERSITY STATISTICS
  E = sum(R, 1) * Pr_b.t();
  O = R * Phi.t();

  compute_objective();
  objective_harmony.push_back(objective_kmeans.back());
  ran_init = true;

}

void harmony::compute_objective() {
  float kmeans_error = as_scalar(accu(R % dist_mat));
  float _entropy = as_scalar(accu(safe_entropy(R).each_col() % sigma)); // NEW: vector sigma
  float _cross_entropy;
  _cross_entropy = as_scalar(accu((R.each_col() % sigma) % ((arma::repmat(theta.t(), K, 1) % log((O + 1) / (E + 1))) * Phi)));
  objective_kmeans.push_back(kmeans_error + _entropy + _cross_entropy);
  objective_kmeans_dist.push_back(kmeans_error);
  objective_kmeans_entropy.push_back(_entropy);
  objective_kmeans_cross.push_back(_cross_entropy);
}


bool harmony::check_convergence(int type) {
  float obj_new, obj_old;
  switch (type) {
  case 0:
    // Clustering
    // compute new window mean
    obj_old = 0;
    obj_new = 0;
    for (int i = 0; i < window_size; i++) {
      obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
      obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
    }
    if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans) {
      return(true);
    } else {
      return(false);
    }
  case 1:
    // Harmony
    obj_old = objective_harmony[objective_harmony.size() - 2];
    obj_new = objective_harmony[objective_harmony.size() - 1];
    if ((obj_old - obj_new) / abs(obj_old) < epsilon_harmony) {
      return(true);
    } else {
      return(false);
    }
  }

  // gives warning if we don't give default return value
  return(true);
}





int harmony::cluster_cpp() {
  int err_status = 0;
  int iter;
  Progress p(max_iter_kmeans, verbose);

  // Z_cos has changed
  // R has assumed to not change
  // so update Y to match new integrated data
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
  for (iter = 0; iter < max_iter_kmeans; iter++) {
    p.increment();
    if (Progress::check_abort())
      return(-1);

    // STEP 1: Update Y
    Y = compute_Y(Z_cos, R);
    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed

    // STEP 3: Update R
    err_status = update_R();
    if (err_status != 0) {
      // Rcout << "Compute R failed. Exiting from clustering." << endl;
      return err_status;
    }

    // STEP 4: Check for convergence
    compute_objective();
    if (iter > window_size) {
      bool convergence_status = check_convergence(0);
      if (convergence_status) {
        //        Rcout << "... Breaking Clustering ..., status = " << convergence_status << endl;
        iter++;
        // Rcout << "Clustered for " << iter << " iterations" << endl;
        break;
      }
    }
  }
  kmeans_rounds.push_back(iter);
  objective_harmony.push_back(objective_kmeans.back());
  return 0;
}




int harmony::update_R() {
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
  _scale_dist = -dist_mat;
  _scale_dist.each_col() /= sigma; // NEW: vector sigma
  _scale_dist.each_row() -= max(_scale_dist, 0);
  _scale_dist = exp(_scale_dist);

  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  int n_blocks = (int)(ceil(1.0 / block_size));
  int cells_per_block = (N / n_blocks) + 1;
  for (int i = 0; i < n_blocks; i++) {
    // gather cell updates indices
    int idx_min = i * cells_per_block;
    int idx_max = min(idx_min + cells_per_block, N);
    if (idx_min > idx_max) break;
    uvec idx_list = linspace<uvec>(idx_min, idx_max - 1, idx_max - idx_min);
    cells_update = update_order.rows(idx_list);

    // Step 1: remove cells
    E -= sum(R.cols(cells_update), 1) * Pr_b.t();
    O -= R.cols(cells_update) * Phi.cols(cells_update).t();

    // Step 2: recompute R for removed cells
    R.cols(cells_update) = _scale_dist.cols(cells_update);
    R.cols(cells_update) = R.cols(cells_update) % (pow((E + 1) / (O + 1), theta) * Phi.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns

    // Step 3: put cells back
    E += sum(R.cols(cells_update), 1) * Pr_b.t();
    O += R.cols(cells_update) * Phi.cols(cells_update).t();

  }
  return 0;
}


void harmony::moe_correct_ridge_cpp() {
  Z_corr = Z_orig;
  for (int k = 0; k < K; k++) {
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
    W.row(0).zeros(); // do not remove the intercept
    Z_corr -= W.t() * Phi_Rk;
  }
  Z_cos = arma::normalise(Z_corr, 2, 0);
}

CUBETYPE harmony::moe_ridge_get_betas_cpp() {
  CUBETYPE W_cube(W.n_rows, W.n_cols, K); // rows, cols, slices
  for (unsigned k = 0; k < K; k++) {
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W_cube.slice(k) = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
  }
  return W_cube;
}

RCPP_MODULE(harmony_module) {
  class_<harmony>("harmony")
  .constructor<int>()

  .field("Z_corr", &harmony::Z_corr)
  .field("Z_orig", &harmony::Z_orig)
  .field("Z_cos", &harmony::Z_cos)
  .field("R", &harmony::R)
  .field("Y", &harmony::Y)
  .field("Phi", &harmony::Phi)
  .field("Phi_moe", &harmony::Phi_moe)
  .field("Pr_b", &harmony::Pr_b)
  .field("objective_kmeans", &harmony::objective_kmeans)
  .field("objective_kmeans_dist", &harmony::objective_kmeans_dist)
  .field("objective_kmeans_entropy", &harmony::objective_kmeans_entropy)
  .field("objective_kmeans_cross", &harmony::objective_kmeans_cross)
  .field("objective_harmony", &harmony::objective_harmony)
  .field("dist_mat", &harmony::dist_mat)
  .field("ran_setup", &harmony::ran_setup)
  .field("ran_init", &harmony::ran_init)

  .field("N", &harmony::N)
  .field("K", &harmony::K)
  .field("B", &harmony::B)
  .field("d", &harmony::d)
  .field("W", &harmony::W)
  .field("max_iter_kmeans", &harmony::max_iter_kmeans)

  .field("sigma", &harmony::sigma)
  .field("theta", &harmony::theta)
  .field("lambda", &harmony::lambda)
  .field("O", &harmony::O)
  .field("E", &harmony::E)
  .field("update_order", &harmony::update_order)
  .field("cells_update", &harmony::cells_update)
  .field("kmeans_rounds", &harmony::kmeans_rounds)
  .field("epsilon_kmeans", &harmony::epsilon_kmeans)
  .field("epsilon_harmony", &harmony::epsilon_harmony)

  // .method("init_cluster", &harmony::init_cluster)
  .method("check_convergence", &harmony::check_convergence)
  .method("setup", &harmony::setup)
  .method("compute_objective", &harmony::compute_objective)
  .method("update_R", &harmony::update_R)
  .method("init_cluster_cpp", &harmony::init_cluster_cpp)
  .method("cluster_cpp", &harmony::cluster_cpp)
  .method("moe_correct_ridge_cpp", &harmony::moe_correct_ridge_cpp)
  .method("moe_ridge_get_betas_cpp", &harmony::moe_ridge_get_betas_cpp)

  ;
}







