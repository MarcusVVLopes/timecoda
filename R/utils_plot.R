# Module 5: helpers shared by the plot methods across modules.
#
# `tc_triplet_plots()` is called from four files, `tc_draw_plots()` from five
# (Modules 0, 1, 2 and 4), so they live here rather than beside any one
# caller. `tc_warn_unrendered()` accompanies `tc_draw_plots()`,
# its only caller. `tc_check_plot_dots()` is called from the two surface
# plots, in two files.

# One ternary panel per triplet of `parts`, built by `panel_fn(tri)` and named
# "a.b.c": the bare plot when there is exactly one triplet, else a named list.
# Shared by every ternary plot method, which is why the three-part minimum is
# enforced here. A plot's coord is set to draw unclipped: ggtern draws a corner
# title inside the panel, so a long one loses its outer letters at the panel's
# edge, and coord_tern() has no `clip` argument to set it with.
tc_triplet_plots <- function(parts, panel_fn, call = rlang::caller_env()) {
  if (length(parts) < 3) {
    cli::cli_abort(
      c(
        "A ternary plot needs at least three composition parts.",
        "i" = "This composition has {length(parts)} part{?s}."
      ),
      call = call
    )
  }
  triplets <- utils::combn(parts, 3, simplify = FALSE)
  plots <- purrr::map(triplets, function(tri) {
    p <- panel_fn(tri)
    if (inherits(p, "gg")) p$coordinates$clip <- "off"
    p
  })
  names(plots) <- purrr::map_chr(triplets, paste, collapse = ".")
  if (length(plots) == 1) plots[[1]] else plots
}

# Draw a ggplot/ggtern object or a (possibly nested) list of them, so that a
# panel failing to render degrades to a warning instead of aborting the rest.
# Returns `x` invisibly, so a plot method can end with this call.
# `ggtern` registers ternary theme elements that some `ggplot2` versions then
# reject when building any plot, so on that error the registered elements are
# cleared (keeping the user's own theme) and the draw is retried once.
tc_draw_plots <- function(x) {
  if (inherits(x, "gg")) {
    tryCatch(
      print(x),
      error = function(e) {
        if (!grepl("tern\\.", conditionMessage(e), perl = TRUE)) {
          return(tc_warn_unrendered(e))
        }
        default_theme <- ggplot2::theme_get()
        ggplot2::reset_theme_settings()
        ggplot2::theme_set(default_theme)
        tryCatch(print(x), error = tc_warn_unrendered)
      }
    )
  } else if (is.list(x)) {
    purrr::walk(x, tc_draw_plots)
  }
  invisible(x)
}

# Warn that a panel could not be drawn, naming the likely ggtern/ggplot2 cause.
tc_warn_unrendered <- function(e) {
  cli::cli_warn(c(
    "Could not render a plot panel.",
    "x" = conditionMessage(e),
    "i" = "If this is a ternary plot, the installed {.pkg ggtern} is likely incompatible with {.pkg ggplot2} {utils::packageVersion('ggplot2')}; update {.pkg ggtern}."
  ), class = "timecoda_unrendered_panel")
}

# Abort when the surface plot `fn` was given an argument it does not take: a
# name falls into `...` and would otherwise be ignored. The arguments are named,
# an unnamed one by its position in `...`; `.hint` is an optional "i" bullet.
# The dotted names keep the caller's stray arguments out of this signature.
tc_check_plot_dots <- function(.fn, ..., .hint = NULL,
                               .call = rlang::caller_env()) {
  n <- ...length()
  if (n == 0) {
    return(invisible())
  }
  nms <- ...names()
  if (is.null(nms)) nms <- rep("", n)
  nms[is.na(nms)] <- ""
  given <- purrr::map2_chr(nms, seq_len(n), function(nm, k) {
    if (nzchar(nm)) {
      cli::format_inline("{.arg {nm}}")
    } else {
      cli::format_inline("an unnamed argument in position {k} of {.code ...}")
    }
  })
  fn <- .fn
  cli::cli_abort(c(
    "{.fn {fn}} was given {cli::qty(n)}{?an argument/arguments} it does not take: {given}.",
    "i" = .hint
  ), call = .call)
}
