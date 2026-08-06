# Column names and the magrittr dot referred to inside NSE calls, declared so
# R CMD check does not read them as undefined globals.
utils::globalVariables(c(".", "x", "y"))
