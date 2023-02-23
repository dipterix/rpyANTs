#' @export
`$.ants.proxy` <- function(x, name) {
  re <- NextMethod("$")
  if(inherits(re, "python.builtin.object")) {
    cls <- class(re)
    if(!"ants.proxy" %in% cls) {
      class(re) <- c("ants.proxy", cls)
    }
  }
  re
}


#' @export
print.ants.proxy <- function(x, ...) {

  str <- reticulate::py_capture_output(
    type = "stdout",
    {
      base <- reticulate::import_builtins()
      base$help(x)
    }
  )
  cat("<ANTs Python Wrapper>\n")
  cat(str)
  cat("\n*** Above is documentation for Python. Please use `$` in R accordingly\n")
  NextMethod()
  invisible(x)
}
