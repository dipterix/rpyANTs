#' @export
`$.ants.proxy` <- function(x, name) {
  re <- NextMethod("$")
  if(inherits(re, c(
    "python.builtin.type",
    "python.builtin.function",
    "python.builtin.module"
  ))) {
    cls <- class(re)
    if(!"ants.proxy" %in% cls) {
      class(re) <- c("ants.proxy", cls)
    }
  }
  re
}


#' @export
print.ants.proxy <- function(x, ...) {

  cat("<ANTs Python Wrapper>\n")
  tryCatch({
    str <- reticulate::py_capture_output(
      type = "stdout",
      {
        base <- reticulate::import_builtins()
        base$help(x)
      }
    )
    str <- trimws(str)
    cat(str)
    cat("\n\n*** Above documentation is for Python. \n*** Please use `$` instead of `.` for modules and functions in R\n")
  }, error = function(e) {
    cat("Cannot retrieve help function...\n")
  })

  NextMethod()
  invisible(x)
}

