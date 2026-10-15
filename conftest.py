# Use a non-interactive matplotlib backend for the test run, so functions that
# call plt.show() (e.g. estimate_read_depth) don't open a window and block on
# a machine with a display. CI is headless and falls back to Agg already.
import matplotlib

matplotlib.use("Agg")
