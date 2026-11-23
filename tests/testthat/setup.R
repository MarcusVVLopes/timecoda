# Several plot functions (plot.compD(), plot.comp_iv(),
# plot_ternary(), plot_realloc()) draw as a side effect. Route that drawing to a
# throwaway device for the whole test run so it does not litter an Rplots.pdf;
# the device is closed automatically when the test session ends.
grDevices::pdf(NULL)
withr::defer(grDevices::dev.off(), teardown_env())
