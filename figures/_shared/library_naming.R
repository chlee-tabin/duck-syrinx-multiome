# Single source of truth for tech-rep library collapse.
#
# Mirrors figures/_shared/library_naming.py::TECH_REP.  The R ggplot
# scripts historically kept their own copies; any tech-rep addition that
# lands in one language but not the other breaks silently.  Source this
# file instead.
#
# Usage:
#   source(file.path(dirname(this_path), "..", "_shared", "library_naming.R"))
#   df$library_m <- collapse_lib(df$library)

TECH_REP <- c(
  "GEX03-1A"   = "GEX03-1",
  "GEX03-1B"   = "GEX03-1",
  "PILOTGEX4A" = "PILOTGEX4",
  "PILOTGEX4B" = "PILOTGEX4"
)

collapse_lib <- function(raw) {
  ifelse(raw %in% names(TECH_REP), unname(TECH_REP[raw]), raw)
}
