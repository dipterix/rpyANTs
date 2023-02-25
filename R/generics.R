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
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__eq__`(e2))
  } else {
    return(e2$`__eq__`(e1))
  }
}

#' @export
`>=.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__ge__`(e2))
  } else {
    return(e2$`__le__`(e1))
  }
}

#' @export
`>.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__gt__`(e2))
  } else {
    return(e2$`__lt__`(e1))
  }
}

#' @export
`<=.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__le__`(e2))
  } else {
    return(e2$`__ge__`(e1))
  }
}

#' @export
`<.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__lt__`(e2))
  } else {
    return(e2$`__gt__`(e1))
  }
}

#' @export
`*.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__mul__`(e2))
  } else {
    return(e2$`__mul__`(e1))
  }
}

#' @export
`!=.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__ne__`(e2))
  } else {
    return(e2$`__ne__`(e1))
  }
}

#' @export
`^.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  e1$`__pow__`(e2)
}

#' @export
`+.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__add__`(e2))
  } else {
    return(e2$`__add__`(e1))
  }
}

#' @export
`-.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__sub__`(e2))
  } else {
    return(e2$`__mul__`(-1L)$`__add__`(e1))
  }
}

#' @export
`/.ants.core.ants_image.ANTsImage` <- function(e1, e2) {
  if(inherits(e1, "ants.core.ants_image.ANTsImage")) {
    return(e1$`__truediv__`(e2))
  } else {
    x <- e2$clone()
    x[] <- e1 / e2[drop = FALSE]
    return(x)
  }
}

#' @export
as.array.ants.core.ants_image.ANTsImage <- function(x, ...) {
  # do NOT use numpy as data is copied anyway from py to r
  py_to_r(x$view())
}

#' @export
dim.ants.core.ants_image.ANTsImage <- function(x) {
  as.integer(py_to_r(x$shape))
}

#' @export
`[.ants.core.ants_image.ANTsImage` <- function(x, i, ..., drop = TRUE) {
  if(!missing(i) && is_py_inherits(i)) {
    x <- py_to_r(x$`__getitem__`(i))
    if(drop) {
      x <- baseenv()$drop(x)
    }
  } else {
    x <- py_to_r(x$view())
    x[i, ..., drop = drop]
  }
}

#' @export
`[<-.ants.core.ants_image.ANTsImage` <- function(x, i, ..., inplace = FALSE, value) {
  if(!inplace) {
    x <- x$clone()
  }
  if(!missing(i) && is_py_inherits(i)) {
    x$`__setitem__`(i, value)
  } else {
    arr <- py_to_r(x$view())
    arr[i, ...] <- value
    x$`__setitem__`(NULL, arr)
  }
  x
}


#' @export
min.ants.core.ants_image.ANTsImage <- function(x, ..., axis = NULL, na.rm = FALSE) {
  axis <- convert_if_not_python(axis, {
    if(length(axis)) {
      as.integer(axis)
    } else {
      NULL
    }
  })
  np <- import("numpy", convert = TRUE)
  np$asarray(x$min(axis))
}

#' @export
max.ants.core.ants_image.ANTsImage <- function(x, ..., axis = NULL, na.rm = FALSE) {
  axis <- convert_if_not_python(axis, {
    if(length(axis)) {
      as.integer(axis)
    } else {
      NULL
    }
  })
  np <- import("numpy", convert = TRUE)
  np$asarray(x$max(axis))
}

#' @export
range.ants.core.ants_image.ANTsImage <- function(x, ..., axis = NULL, na.rm = FALSE) {
  axis <- convert_if_not_python(axis, {
    if(length(axis)) {
      as.integer(axis)
    } else {
      NULL
    }
  })
  np <- import("numpy", convert = TRUE)
  np$asarray(x$range(axis))
}

#' @export
log.ants.core.ants_image.ANTsImage <- function(x, base = exp(1)) {

  y <- log(py_to_r(x$view()), base = base)
  x <- x$clone()
  x[] <- y
  x

}

#' @export
exp.ants.core.ants_image.ANTsImage <- function(x) {

  y <- exp(py_to_r(x$view()))
  x <- x$clone()
  x[] <- y
  x

}
