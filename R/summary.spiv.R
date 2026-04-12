#' @title Summary for SP-IV models
#'
#' @description Extracts a summary for \code{\link{spiv}} class models.
#'
#' @param object A `spiv` object.
#'
#' @return Prints a formatted text summary of the estimates and diagnostics.
#'
#' @export

summary.spiv <- function(object, ...) {
  cat("\n=== SP-IV Estimator Results ===\n")

  cat("\n--- Structural Estimates (Strong IV) ---\n")
  print(object$strong_ci[, c("Estimate", "Std_Error", "Lower_95", "Upper_95")], digits = 4)

  cat("\n--- Weak Instrument Diagnostics ---\n")
  cat("Statistic:      ", round(object$weak_iv_test$statistic, 3), "\n")
  cat("Critical Value: ", round(object$weak_iv_test$critical_value, 3), "\n")
  cat("Instruments Weak?", object$weak_iv_test$is_weak, "\n")

  cat("\n--- Robust Inference (", object$robust_inference$method, ") ---\n", sep="")
  if (nrow(object$robust_inference$bounds) > 0) {
    print(object$robust_inference$bounds, digits = 4)
  } else {
    cat("WARNING:", object$robust_inference$warning, "\n")
  }
  cat("===============================\n\n")
}
