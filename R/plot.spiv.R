#' @title Plot for SP-IV models
#' 
#' @description Draws a plot for \code{\link{spiv}} class models.
#' 
#' @param x A `spiv` object.
#' @param type The IRF to plot. `"outcome"` gives the IRF of the outcome variable $y$ and `"endogenous"` gives the IRFs of the endogenous variables.
#' @param idx If `type = c"outcome"`, `idx` must be the index of the instrument you want the IRF with respect to or `all` if you want a grid with every plot. If `type = "endogenous"`, `idx` must be a column vector : `c(endogenous variable index, instrument index)` or `all`.
#' @param ... Additional arguments. 
#' 
#' @return A ggplot of the IRF and HAC confidence intervals.
#' 
#' @importFrom ggplot2 ggplot aes geom_hline geom_ribbon geom_line theme_minimal labs facet_wrap facet_grid
#' @importFrom tools toTitleCase
#'
#' @export

plot.spiv <- function(x, type = c("outcome", "endogenous"), idx = "all", ...) {
  object <- x
  type <- match.arg(type)
  dat <- object$irf[[type]]
  
  # Extract dimensions from the object
  H <- object$args$H
  Nz <- object$generated_args$Nz
  K <- object$generated_args$K
  
  # Grid plot
  if (length(idx) == 1 && idx == "all") {
    
    df_list <- list()
    counter <- 1
    
    if (type == "outcome") {
      for (s in 1:Nz) {
        df_list[[counter]] <- data.frame(
          Horizon = 0:(H - 1),
          Estimate = dat$point_estimate[, s],
          Lower = dat$lower_bound[, s],
          Upper = dat$upper_bound[, s],
          Variable = "Outcome",
          Shock = paste("Shock", s)
        )
        counter <- counter + 1
      }
    } else {
      for (v in 1:K) {
        start_row <- (v - 1) * H + 1
        end_row <- v * H
        for (s in 1:Nz) {
          df_list[[counter]] <- data.frame(
            Horizon = 0:(H - 1),
            Estimate = dat$point_estimate[start_row:end_row, s],
            Lower = dat$lower_bound[start_row:end_row, s],
            Upper = dat$upper_bound[start_row:end_row, s],
            Variable = paste("Endogenous", v),
            Shock = paste("Shock", s)
          )
          counter <- counter + 1
        }
      }
    }
    
    # Combine the list into one dataframe for ggplot
    df_all <- do.call(rbind, df_list)
    
    # Plot with facets
    g <- ggplot(df_all, aes(x = Horizon, y = Estimate)) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "gold3") +
      geom_line(color = "gold3", linewidth = 1) +
      theme_minimal() +
      labs(
        title = paste("SP-IV Impulse Responses:", tools::toTitleCase(type)),
        x = "Horizon (h)",
        y = "Response"
      )
    
    if (type == "outcome") {
      # 1 Row, Nz Columns
      g <- g + facet_wrap(~ Shock, ncol = Nz, scales = "free_y")
    } else {
      # K Rows, Nz Columns (Automatically forms the exact Theta_Y matrix shape)
      g <- g + facet_grid(Variable ~ Shock, scales = "free_y")
    }
    
    return(g)
  }
  
  # Single plot 
  if (type == "outcome") {
    shock_idx <- idx[1]
    if (shock_idx > Nz) stop(paste("ERROR: Shock index exceeds number of instruments (", Nz, ")."))
    
    est <- dat$point_estimate[, shock_idx]
    lw  <- dat$lower_bound[, shock_idx]
    up  <- dat$upper_bound[, shock_idx]
    
    plot_title <- "Impulse Response: Outcome"
    plot_sub <- paste("Response to Shock", shock_idx)
    
  } else if (type == "endogenous") {
    if (length(idx) < 2) stop("ERROR: For specific 'endogenous' plots, 'idx' must be c(variable_index, shock_index). Or use idx='all'.")
    
    var_idx <- idx[1]
    shock_idx <- idx[2]
    
    if (var_idx > K) stop(paste("ERROR: Variable index exceeds total endogenous variables (", K, ")."))
    if (shock_idx > Nz) stop(paste("ERROR: Shock index exceeds total instruments (", Nz, ")."))
    
    start_row <- (var_idx - 1) * H + 1
    end_row <- var_idx * H
    
    est <- dat$point_estimate[start_row:end_row, shock_idx]
    lw  <- dat$lower_bound[start_row:end_row, shock_idx]
    up  <- dat$upper_bound[start_row:end_row, shock_idx]
    
    plot_title <- paste("SP-IV Impulse Response: Endogenous Variable", var_idx)
    plot_sub <- paste("Response to Shock", shock_idx)
  }
  
  df <- data.frame(Horizon = 0:(H - 1), Estimate = est, Lower = lw, Upper = up)
  
  ggplot(df, aes(x = Horizon, y = Estimate)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "pink4") +
    geom_line(color = "pink4", linewidth = 1) +
    theme_minimal() +
    labs(
      title = plot_title,
      subtitle = plot_sub,
      x = "Horizon (h)",
      y = "Response"
    )
}