# CytochromeP450-specific constants and functions

HELIX_ORDER = c("A'", "A", "B", "B'", "B''", "C", "D", "E", "F", "F'", "G'", "G", "H", "I", "J", "J'", "K", "K'", "K''", "L", "L'")
MAJOR_HELICES = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
SHEET_ORDER = c("1-0", "1-1", "1-2", "1-3", "1-4", "1-5", "2-1", "2-2", "3-1", "3-2", "3-3", "4-1", "4-2", "5-1", "5-2", "6-1", "6-2")
OUR_SSES = c(HELIX_ORDER, SHEET_ORDER)

# Decide if the SSE label corresponds to major, minor helix, or strand
sse_category = function(name) {
  name = as.character(name)
  if (name %in% MAJOR_HELICES) 'major helix'
  else if (name %in% HELIX_ORDER) 'minor helix'
  else if (name %in% SHEET_ORDER) 'strand'
  else 'other'
}
