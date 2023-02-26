
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
  to_r(x$view())
}

#' @export
dim.ants.core.ants_image.ANTsImage <- function(x) {
  as.integer(to_r(x$shape))
}

#' @export
`[.ants.core.ants_image.ANTsImage` <- function(x, i, ..., drop = TRUE) {
  if(!missing(i) && is_py_inherits(i)) {
    x <- to_r(x$`__getitem__`(i))
    if(drop) {
      x <- baseenv()$drop(x)
    }
  } else {
    x <- to_r(x$view())
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
    arr <- to_r(x$view())
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

  y <- log(to_r(x$view()), base = base)
  x <- x$clone()
  x[] <- y
  x

}

#' @export
exp.ants.core.ants_image.ANTsImage <- function(x) {

  y <- exp(to_r(x$view()))
  x <- x$clone()
  x[] <- y
  x

}
