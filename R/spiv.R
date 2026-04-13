#' @title System Projections with Instrumental Variables
#'
#' @description Estimates structural parameters and dynamic impulse responses using the SP-IV methodology (D. J. Lewis and K. Mertens, 2024, available at \url{https://karelmertens.com/wp-content/uploads/2024/07/dynamic_identification_using_system_projections_and_instrumental_variables_main_7_1_24.pdf}). The paper used to source the equations is the Federal Reserve Bank of Dallas working paper from March 2022 (\doi{https://doi.org/10.24149/wp2204}). The implementation uses Local Projections (section 4.1 of the Dallas Fed WP). Both AR and KLM tests are available but KLM can be very slow for large grids.
#'
#' @param y A T x 1 matrix representing the outcome variable.
#' @param Y A T x K matrix representing the endogenous variables.
#' @param X A Nx x T matrix representing the control variables.
#' @param Z A Nz x T matrix representing the instrumental variables.
#' @param H An integer specifying the maximum forecast horizon.
#' @param alpha A numeric value for the significance level (default is `0.05`).
#' @param wiv A character string specifying the weak instrument test to use: `"AR"` for Anderson-Rubin or `"KLM"` for Kleibergen.
#' @param xi A numeric value specifying the bias tolerance for the weak instrument test (default is `0.10`) : see section 4.3.2 of the Dallas Fed WP.
#' @param grid A named numeric vector (lower, upper, length) defining the bounds and resolution of the robust grid search (default is `c(-1, 1, 100)`).
#' @param cores An integer specifying the number of CPU cores to use for parallel processing during the grid search (default is `1`).
#'
#' @return An S3 object of class \code{spiv} containing point estimates, standard errors, weak instrument diagnostics, robust bounds, impulse response matrices and arguments used.
#'
#' @useDynLib spiv, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom sandwich NeweyWest
#' @importFrom lmtest coeftest
#' @importFrom stats qnorm qchisq pnorm
#'
#' @examples
#' \dontrun{
#' # Assuming y, Y, X, and Z are pre-loaded matrices
#' results <- spiv(y = y, Y = Y, X = X, Z = Z, H = 10, wiv = "AR", cores = 4)
#' summary(results)
#' plot(results)
#' }
#'
#' @export
spiv <- function(y, Y, X, Z, H, alpha = 0.05, wiv = c("AR", "KLM"), xi = 0.10, grid = c(lower = -1, upper = 1, length = 100), cores = 1) {
    plan(multisession, workers = cores)
    wiv <- match.arg(wiv)
    inputs_list <- mget(names(formals()))

    if (!is.matrix(y) || !is.matrix(Y) || !is.matrix(X) || !is.matrix(Z)) {
      stop("ERROR: Input data (y, Y, X, Z) must be matrices. Use as.matrix() if passing data frames.")
    }

    T_total <- nrow(Y)
    if (nrow(y) != T_total) {
      stop("ERROR: Dimension mismatch. 'y' and 'Y' must have the same number of rows.")
    }
    if (ncol(X) != T_total) {
      stop("ERROR: Dimension mismatch. 'X' must have the same number of columns as the rows of 'Y'.")
    }
    if (ncol(Z) != T_total) {
      stop("ERROR: Dimension mismatch. 'Z' must have the same number of columns as the rows of 'Y'.")
    }

    if (H < 1 || H >= T_total) {
      stop("ERROR: Horizon 'H' must be a positive integer strictly less than the total time periods.")
    }

    if (H * nrow(Z) < ncol(Y)) {
      stop("ERROR: Under-identified. The SP-IV order condition requires H * Nz >= K.")
    }

  # Pre-process
    K <- ncol(Y)
    Nx <- nrow(X)
    Nz <- nrow(Z)
    T_ess <- nrow(Y) - H + 1 # Sample size accounting for leads

    if (T_ess - Nx - K <= 0) {
      stop("ERROR: Insufficient degrees of freedom for the structural error variance (T_ess - Nx - K <= 0).")
    }
    if (T_ess - Nx - Nz <= 0) {
      stop("ERROR: Insufficient degrees of freedom for the first-stage error variance (T_ess - Nx - Nz <= 0).")
    }

    Y_hlist <- lapply(1:K, function(k) {
      leads <- sapply(0:(H-1), function(h) {
        Y[(1+h):(T_ess + h), k] # This is Y_k^H
      })
      return(t(leads))
    })
    YH <- do.call(rbind, Y_hlist)

    y_hlist <- sapply(0:(H-1), function(h) {
      y[(1+h):(T_ess + h), 1]
    })
    yH <- t(y_hlist)

    # Truncate X and Z to match the new T_ess dimension
    X <- X[, 1:T_ess]
    Z <- Z[, 1:T_ess]

  # Direct forecasting
    Px = t(X) %*% qr.solve(X %*% t(X)) %*% X
    Mx = diag(T_ess) - Px

    # Equation 38 :
    y_res = yH %*% Mx
    Y_res = YH %*% Mx
    Z_res = Z %*% Mx

  # SP-IV estimator
    R <- diag(K) %x% as.vector(diag(H))
    Pz <- t(Z_res) %*% solve(Z_res %*% t(Z_res)) %*% Z_res
    Mz <- diag(T_ess) - Pz

    # Eq 9 :
    beta <- qr.solve(t(R) %*% ((Y_res %*% Pz %*% t(Y_res)) %x% diag(H)) %*% R) %*% t(R) %*% as.vector(y_res %*% Pz %*% t(Y_res))

  # As an IRF (section 2.4)
    eig <- eigen(qr.solve((Z_res %*% t(Z_res)) / T_ess))
    inv_sqrt_Z <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

    Theta_Y <- ((Y_res %*% t(Z_res))/T_ess) %*% inv_sqrt_Z
    Theta_y <- ((y_res %*% t(Z_res))/T_ess) %*% inv_sqrt_Z

  # Error variance (4.3.1)
    u_res <- y_res - (t(beta) %x% diag(H)) %*% Y_res

    v_eig <- eigen(qr.solve(Z %*% Mx %*% t(Z) / T_ess))
    v_inv_sqrtm <- v_eig$vectors %*% diag(sqrt(v_eig$values)) %*% t(v_eig$vectors)
    v_res <- Y_res - Theta_Y %*% v_inv_sqrtm %*% Z %*% Mx

    vcov_uH <- (u_res %*% t(u_res))/(T_ess-Nx-K)
    vcov_vH <- (v_res %*% t(v_res))/(T_ess-Nx-Nz)

    vcov_vH_eig <- eigen(vcov_vH)
    vcov_vH_sqrtm <- vcov_vH_eig$vectors %*% diag(sqrt(vcov_vH_eig$values)) %*% t(vcov_vH_eig$vectors)

  # Testing for Weak Instruments
    # 1
    Omega <- t(R) %*% (vcov_vH %x% diag(H)) %*% R
    Omega_eig <- eigen(qr.solve(Omega))
    Omega_inv_sqrtm <- Omega_eig$vectors %*% diag(sqrt(Omega_eig$values)) %*% t(Omega_eig$vectors)

    S <- (Omega_inv_sqrtm %x% diag(H)) %*% vcov_vH_sqrtm
    Scal <- S %*% t(S)

    # 2
    g_matrix <- Omega_inv_sqrtm %*% t(R) %*% (Y_res %*% Pz %*% t(Y_res) %x% diag(H)) %*% R %*% Omega_inv_sqrtm
    g_eigen <- eigen(g_matrix)
    g_min <- min(g_eigen$values)/Nz

    # 3
    # function arguments
    ell <- 1/xi

    maxeval_square_eigen <- eigen(t(R) %*% ((Scal %*% Scal) %x% diag(H)) %*% R)
    maxeval_square <- max(maxeval_square_eigen$values)
    maxeval_cube_eigen <- eigen(t(R) %*% ((Scal %*% Scal %*% Scal) %x% diag(H)) %*% R)
    maxeval_cube <- max(maxeval_cube_eigen$values)

    kappa_1 <- Nz * (1+ell)
    kappa_2 <- 2 * (Nz * maxeval_square + 2 * ell * Nz)
    kappa_3 <- 8 * (Nz * maxeval_cube + 3 * ell * Nz * maxeval_square)

    nu <- kappa_2 / kappa_3
    delta <- 8 * kappa_2 * nu^2
    g_crit <- ((qchisq(1-alpha, delta) - delta)/(4 * nu) + kappa_1) / Nz

  # Strong instrument confidence sets (4.3.3)
    V_beta <- (qr.solve(t(R) %*% (Theta_Y %*% t(Theta_Y) %x% diag(H)) %*% R) %*% t(R) %*% (Theta_Y %*% t(Theta_Y) %x% vcov_uH) %*% R %*% qr.solve(t(R) %*% (Theta_Y %*% t(Theta_Y) %x% diag(H)) %*% R))/T_ess

    se_beta <- sqrt(diag(V_beta))

    strong_lower <- beta - qnorm(1-alpha/2) * se_beta
    strong_upper <- beta + qnorm(1-alpha/2) * se_beta

    # output
    strong_ci <- data.frame(
      Estimate = as.vector(beta),
      Std_Error = as.vector(se_beta),
      Lower_95 = as.vector(strong_lower),
      Upper_95 = as.vector(strong_upper)
    )
    rownames(strong_ci) <- paste0("Param_", 1:K)

  # Robust tests
    if (wiv[1] == "AR") {
      # AR
      d_AR = Nz + Nx
      crit_AR <- qchisq(1-alpha, df = H * Nz)

      # Dynamic grid for K variables : change the from, to and length to fct arguments later
      grid_list <- replicate(K, seq(from = grid["lower"], to = grid["upper"], length.out = grid["length"]), simplify = FALSE)
      # Generate the combinations dynamically using do.call
      b_grid_df <- do.call(expand.grid, grid_list)
      # Data frame to list of transposed Kx1 matrices
      b_list <- lapply(1:nrow(b_grid_df), function(i) t(as.matrix(b_grid_df[i, ])))
      # Grid search
      results <- future_lapply(b_list, function(b_current) {
        stat <- compute_AR_stat(b_current, y_res, Y_res, Pz, Mz, H, T_ess, d_AR)
        if (stat < crit_AR) return(b_current) else return(NULL)
      }, future.seed = TRUE)
      # Filter the NULL
      confidence_set <- Filter(Negate(is.null), results)
      # data frame
      if (length(confidence_set) > 0) {
        valid_b_df <- do.call(rbind, lapply(confidence_set, t))
        colnames(valid_b_df) <- paste0("Param_", 1:K)
      } else {
        valid_b_df <- data.frame()
      }
    } else if (wiv[1] == "KLM") {
      # KLM
      d_K = Nz + Nx
      crit_KLM <- qchisq(1-alpha, df = K)

      # Dynamic grid for K variables : change the from, to and length to fct arguments later
      grid_list <- replicate(K, seq(from = grid["lower"], to = grid["upper"], length.out = grid["length"]), simplify = FALSE)
      # Generate the combinations dynamically using do.call
      b_grid_df <- do.call(expand.grid, grid_list)
      # Data frame to list of transposed Kx1 matrices
      b_list <- lapply(1:nrow(b_grid_df), function(i) t(as.matrix(b_grid_df[i, ])))

      # Grid search with RcppArmadillo
      results <- future_lapply(b_list, function(b_current) {
        stat <- compute_KLM_stat(b_current, y_res, Y_res, v_res, Mz, Pz, R, H, T_ess, d_K)
        if (stat < crit_KLM) return(b_current) else return(NULL)
      }, future.seed = TRUE)
      # Filter the NULL
      confidence_set <- Filter(Negate(is.null), results)
      # data frame
      if (length(confidence_set) > 0) {
        valid_b_df <- do.call(rbind, lapply(confidence_set, t))
        colnames(valid_b_df) <- paste0("Param_", 1:K)
      } else {
        valid_b_df <- data.frame()
      }

    }
    # Empty grid check
    if (nrow(valid_b_df) > 0) {
      # inference results
      robust_inference = list(
        method = wiv,
        bounds = data.frame(
          Parameter = colnames(valid_b_df),
          Lower = apply(valid_b_df, 2, min),
          Upper = apply(valid_b_df, 2, max),
          row.names = NULL
        ),
        confidence_set = valid_b_df
      )

      if (wiv == "AR") {
        robust_inference$degrees_freedom = d_AR
        robust_inference$critical_value = crit_AR
      } else if (wiv == "KLM") {
        robust_inference$degrees_freedom = d_K
        robust_inference$critical_value = crit_KLM
      }
    } else {
      robust_inference = list(
        warning = "Empty robust confidence set. Reexamine model specification or choice of instruments.",
        method = wiv,
        bounds = data.frame(),
        confidence_set = data.frame()
      )
    }

  # IRF : HAC se
    Z_std <- inv_sqrt_Z %*% Z_res
    # HAC for each horizon
      # Outcome
      se_Theta_y <- matrix(NA, nrow = H, ncol = Nz)
      for (h in 1:H) {
        y_h <- y_res[h,]
        lm_h <- lm(y_h ~ t(Z_std) - 1) # no intercept because already residualized
        hac_results <- lmtest::coeftest(lm_h, vcov = sandwich::NeweyWest(lm_h, lag = h - 1, prewhite = FALSE))
        se_Theta_y[h, ] <- hac_results[, 2] # se's
      }
      # Endogenous
      se_Theta_Y <- matrix(NA, nrow = H * K, ncol = Nz)
      for (h in 1:(H*K)) {
        Y_h <- Y_res[h,]
        lm_h <- lm(Y_h ~ t(Z_std) - 1)
        idx_h <- (h - 1) %% H # Correct horizon for each variable since YH is stacked
        hac_results <- lmtest::coeftest(lm_h, vcov = sandwich::NeweyWest(lm_h, lag = idx_h, prewhite = FALSE))
        se_Theta_Y[h, ] <- hac_results[, 2] # se's
      }
      # Output
      irf = list(
        outcome = list(
          point_estimate = Theta_y,
          std_error = se_Theta_y,
          lower_bound = Theta_y - qnorm(1-alpha/2) * se_Theta_y,
          upper_bound = Theta_y + qnorm(1-alpha/2) * se_Theta_y,
          p = 2 * (1 - pnorm(abs(Theta_y / se_Theta_y)))
          ),
        endogenous = list(
          point_estimate = Theta_Y,
          std_error = se_Theta_Y,
          lower_bound = Theta_Y - qnorm(1-alpha/2) * se_Theta_Y,
          upper_bound = Theta_Y + qnorm(1-alpha/2) * se_Theta_Y,
          p = 2 * (1 - pnorm(abs(Theta_Y / se_Theta_Y)))
          )
      )


    # output
    weak_iv_test = list(
      statistic = g_min,
      critical_value = g_crit,
      is_weak = g_min < g_crit, # true if weak
      weak_stats = list(
        ell = ell,
        kappa_1 = kappa_1,
        kappa_2 = kappa_2,
        kappa_3 = kappa_3,
        nu = nu,
        delta = delta
      )
    )

    generated_args = list(
      T_ess = T_ess,
      K = K,
      Nx = Nx,
      Nz = Nz,
      y_res = y_res,
      Y_res = Y_res,
      Z_res = Z_res,
      R = R,
      vcov_uH = vcov_uH,
      vcov_vH = vcov_vH
    )

    out <- list(
      beta = as.vector(beta),
      vcov = V_beta,
      strong_ci = strong_ci,
      weak_iv_test = weak_iv_test,
      robust_inference = robust_inference,
      irf = irf,
      args = inputs_list,
      generated_args = generated_args
    )

    class(out) <- "spiv"
    return(out)
}
