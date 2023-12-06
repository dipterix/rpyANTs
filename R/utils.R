#' @title Get 'Python' built-in object
#' @param name object name
#' @param convert see \code{\link[reticulate]{import_builtins}}
#' @returns A python built-in object specified by \code{name}
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'
#' # ------ Basic case: use python `int` as an R function ---------
#' py_int <- py_builtin("int")
#'
#' # a is an R object now
#' a <- py_int(9)
#' print(a)
#' class(a)
#'
#' # ------ Use python `int` as a Python function -----------------
#' py_int2 <- py_builtin("int", convert = FALSE)
#'
#' # b in a python object
#' b <- py_int2(9)
#'
#' # There is no '[1] ' when printing
#' print(b)
#' class(b)
#'
#' # convert to R object
#' py_to_r(b)
#'
#'
#'
#' }
#'
#' @export
py_builtin <- function(name, convert = TRUE) {
  builtins <- reticulate::import_builtins(convert = convert)
  builtins[[name]]
}

#' @title Slice index in 'Python' arrays
#' @param ... passing to \code{slice} ('Python')
#' @returns Index slice instance
#'
#' @examples
#'
#'
#' if(interactive() && ants_available()) {
#'
#'   x <- np_array(array(seq(20), c(4, 5)))
#'
#'   # equivalent to x[::2]
#'   x[py_slice(NULL, NULL, 2L)]
#'
#' }
#' @export
py_slice <- function(...) {
  reticulate::import_builtins(convert = FALSE)$slice(...)
}

#' @title List in 'Python'
#' @param ... passing to \code{list} ('Python')
#' @param convert whether to convert the results back into R; default is no
#' @returns List instance, or an R vector if converted
#'
#' @examples
#'
#'
#' if(interactive() && ants_available()) {
#'
#'   py_list(list(1,2,3))
#'   py_list(c(1,2,3))
#'
#'   py_list(array(1:9, c(3,3)))
#'   py_list(list(list(1:3), letters[1:3]))
#'
#' }
#' @export
py_list <- function(..., convert = FALSE) {
  reticulate::import_builtins(convert = convert)$list(...)
}


#' Load data as \code{'ANTsImage'} class
#' @param x data to be converted; this can be an \code{'ANTsImage'} instance,
#' character, \code{'oro.nifti'} object, \code{'niftiImage'} from
#' package \code{'RNifti'}, or \code{'threeBrain.nii'} from package
#' \code{'threeBrain'}
#' @param strict whether \code{x} should not be \code{NULL}
#' @returns An \code{'ANTsImage'} instance; use \code{ants$ANTsImage} to see
#' the 'Python' documentation
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   ants <- load_ants()
#'
#'   # Python string
#'   x1 <- ants$get_ants_data('r16')
#'   as_ANTsImage( x1 )
#'
#'   # R character
#'   nii_path <- system.file(package = "RNifti",
#'                           "extdata", "example.nii.gz")
#'   as_ANTsImage( nii_path )
#'
#'   # niftiImage object
#'   x2 <- RNifti::readNifti(nii_path)
#'   as_ANTsImage( x2 )
#'
#' }
#'
#' @export
as_ANTsImage <- function(x, strict = FALSE) {
  UseMethod("as_ANTsImage")
}

#' @export
as_ANTsImage.ants.core.ants_image.ANTsImage <- function(x, strict = TRUE) {
  return(x)
}

#' @export
as_ANTsImage.oro.nifti <- function(x, strict = TRUE) {
  tfile <- tempfile(fileext = ".nii.gz", pattern = "as_ANTsImage_")
  RNifti::writeNifti(x, file = tfile)
  on.exit({ unlink(tfile) })
  ants <- load_ants()
  return(ants$image_read(tfile))
}

#' @export
as_ANTsImage.niftiImage <- as_ANTsImage.oro.nifti

#' @export
as_ANTsImage.threeBrain.nii <- function(x, strict = TRUE) {
  as_ANTsImage.oro.nifti(x$header)
}

#' @export
as_ANTsImage.character <- function(x, strict = TRUE) {
  if(length(x) != 1 || is.na(x) || trimws(x) == "") {
    if(strict) {
      stop("as_ANTsImage: for string `x`, length(x) must equal to 1 and cannot be empty/NA under strict mode")
    }
    return(NULL)
  }
  ants <- load_ants()
  return(ants$image_read(x))
}

#' @export
as_ANTsImage.python.builtin.str <- function(x, strict = TRUE) {
  as_ANTsImage.character(py_to_r(x), strict = strict)
}

#' @export
as_ANTsImage.default <- function(x, strict = TRUE) {
  if(is.null(x)) {
    if( strict ) {
      stop("as_ANTsImage: input x (image) cannot be NULL under strict mode")
    }
    return(x)
  }
  stop("as_ANTsImage: unsupported input format: [", paste0(class(x), collapse = ","), "].")
}


as_hexcolor <- function(x, ..., use_alpha = FALSE) {
  s <- grDevices::adjustcolor(col = x, ...)
  if(!use_alpha) {
    s <- substr(s, 1, 7)
  }
  s
}


convert_if_not_python <- function(x, value, convert = TRUE) {
  value <- substitute(value)
  if(!inherits(x, "python.builtin.object")) {
    parent_frame <- parent.frame()
    x <- eval(value, parent_frame)
  }
  if(convert && !inherits(x, "python.builtin.object")) {
    x <- r_to_py(x)
  }
  return(x)
}

is_py_inherits <- function(x, class = NULL) {
  inherits(x, c("python.builtin.object", class))
}


to_r <- function(x) {
  tryCatch({
    reticulate::py_to_r(x)
  }, error = function(e) {
    x
  })
}

snapshot_tempfiles <- function() {
  if(!ants_available()) { return() }
  ants <- load_ants()
  injected <- py_builtin("hasattr", convert = TRUE)(ants$utils, "rpyANTsInjected")
  if(!injected) { return() }
  tryCatch({
    to_r(ants$utils$rpyANTsInjection$requested_tempfiles)
  }, error = function(e) {
    NULL
  })
}

remove_tmpfiles <- function(x, ...) {

  tfiles <- py_list(convert = FALSE)
  if(ants_available()) {
    ants <- load_ants()

    injected <- py_builtin("hasattr", convert = TRUE)(ants$utils, "rpyANTsInjected")
    if(injected && isTRUE(to_r(ants$utils$rpyANTsInjected))) {
      tfiles <- ants$utils$rpyANTsInjection$requested_tempfiles
    }
  }

  for(f in x) {
    if(file.exists(f)) {
      unlink(f, ...)
      tryCatch({
        idx <- to_r(tfiles$index(f))
        if(idx > -1) {
          tfiles$remove(f)
        }
      }, error = function(e) {

      })
    }
  }
}


get_os <- function(){
  if("windows" %in% tolower(.Platform$OS.type)){
    return("windows")
  }
  os <- tolower(R.version$os)
  if(startsWith(os, "darwin")){
    return('darwin')
  }
  if(startsWith(os, "linux")){
    return('linux')
  }
  if(startsWith(os, "solaris")){
    return('solaris')
  }
  if(startsWith(os, "win")){
    return('windows')
  }
  return('unknown')
}


normalize_path <- function(path, must_work = NA) {
  path <- unlist(lapply(path, function(p) {
    if(!file.exists(p)) {
      dname <- dirname(p)
      dname <- normalizePath(dname, winslash = "/", mustWork = must_work)
      p <- file.path(dname, basename(p), fsep = "/")
    } else {
      p <- normalizePath(p, winslash = "/", mustWork = must_work)
    }
    p
  }))

  gsub("[/|\\\\]+", "/", path)
}

file_path <- function(...) {
  path <- file.path(..., fsep = "/")
  gsub("[/|\\\\]+", "/", path)
}

dir_create2 <- function(x, showWarnings = FALSE, recursive = TRUE, check = TRUE, ...) {
  if (!dir.exists(x)) {
    dir.create(x, showWarnings = showWarnings, recursive = recursive, ...)
  }
  if (check && !dir.exists(x)) {
    stop('Cannot create directory at ', shQuote(x))
  }
  invisible(normalize_path(x))
}


R_user_dir <- function (package, which = c("data", "config", "cache")) {
  stopifnot(is.character(package), length(package) == 1L)
  which <- match.arg(which)
  home <- normalizePath("~")
  path <- switch(
    which,
    data = {
      if (nzchar(p <-
                 Sys.getenv("R_USER_DATA_DIR")))
        p
      else if (nzchar(p <-
                      Sys.getenv("XDG_DATA_HOME")))
        p
      else if (.Platform$OS.type ==
               "windows")
        file.path(Sys.getenv("APPDATA"), "R",
                  "data")
      else if (Sys.info()["sysname"] == "Darwin")
        file.path(home,
                  "Library", "Application Support", "org.R-project.R")
      else
        file.path(home,
                  ".local", "share")
    },
    config = {
      if (nzchar(p <-
                 Sys.getenv("R_USER_CONFIG_DIR")))
        p
      else if (nzchar(p <-
                      Sys.getenv("XDG_CONFIG_HOME")))
        p
      else if (.Platform$OS.type ==
               "windows")
        file.path(Sys.getenv("APPDATA"), "R",
                  "config")
      else if (Sys.info()["sysname"] == "Darwin")
        file.path(home,
                  "Library", "Preferences", "org.R-project.R")
      else
        file.path(home,
                  ".config")
    },
    cache = {
      if (nzchar(p <-
                 Sys.getenv("R_USER_CACHE_DIR")))
        p
      else if (nzchar(p <-
                      Sys.getenv("XDG_CACHE_HOME")))
        p
      else if (.Platform$OS.type ==
               "windows")
        file.path(Sys.getenv("LOCALAPPDATA"),
                  "R", "cache")
      else if (Sys.info()["sysname"] == "Darwin")
        file.path(home,
                  "Library", "Caches", "org.R-project.R")
      else
        file.path(home,
                  ".cache")
    }
  )
  file.path(path, "R", package)
}
