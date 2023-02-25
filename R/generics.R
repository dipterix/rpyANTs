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


#' @export
`==.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__eq__`(e2)
}

#' @export
`>=.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__ge__`(e2)
}

#' @export
`>.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__gt__`(e2)
}

#' @export
`<=.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__le__`(e2)
}

#' @export
`<.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__lt__`(e2)
}

#' @export
`*.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__mul__`(e2)
}

#' @export
`!=.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__ne__`(e2)
}

#' @export
`^.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__pow__`(e2)
}

#' @export
`+.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__add__`(e2)
}

#' @export
`-.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__sub__`(e2)
}

#' @export
`/.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__truediv__`(e2)
}

