# Draw-time warnings of the package's own class. On a machine where ggtern
# cannot render against the installed ggplot2, every ternary panel warns with
# class "timecoda_unrendered_panel". These tests assert the returned objects,
# so that one class is silenced; any other warning still surfaces.
quiet_draw <- function(expr) {
  suppressWarnings(expr, classes = "timecoda_unrendered_panel")
}

# DHARMa's residual plot fits quantile smooths (qgam, mgcv) whose optimizer
# stops early on some platforms and says so in a warning; the plot is drawn
# all the same. That one warning is silenced; any other still surfaces.
quiet_smooth <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("step failure", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# A d-part simulated composition for the plot tests: 80 rows of p1..pd,
# uniform on [5, 50], and an outcome y, drawn from `seed`.
sim_parts <- function(d, seed = 11) {
  set.seed(seed)
  sim <- as.data.frame(matrix(stats::runif(80 * d, 5, 50), ncol = d))
  names(sim) <- paste0("p", seq_len(d))
  sim$y <- stats::rnorm(80)
  sim
}
