# CytochromeP450-specific constants and functions

HELIX_ORDER = c("A'", "A", "B", "B'", "B''", "C", "D", "E", "F", "F'", "G'", "G", "H", "I", "J", "J'", "K", "K'", "K''", "L", "L'")
SHEET_ORDER = c("1-0", "1-1", "1-2", "1-3", "1-4", "1-5", "2-1", "2-2", "3-1", "3-2", "3-3", "4-1", "4-2", "5-1", "5-2", "6-1", "6-2")
OUR_SSES = c(HELIX_ORDER, SHEET_ORDER)

# Decide if the SSE label corresponds to major, minor helix, or strand
sse_category = function(name) {
  name = as.character(name)
  n = nchar(name)
  if (is_digit(substr(name, 1, 1))) 'strand'
  else if (substr(name, n, n)=='\'') 'minor h.' 
  else 'major h.'
}
