#' @title Get 'Python' built-in object
#' @param name object name
#' @param convert see \code{\link[reticulate]{import_builtins}}
#' @returns A python built-in object
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
#' @export
as_ANTsImage <- function(x, strict = FALSE) {
  if(is.null(x)) {
    if( strict ) {
      stop("as_ANTsImage: input x (image) cannot be NULL under strict mode")
    }
    return(x)
  }
  if(inherits(x, "ants.core.ants_image.ANTsImage")) {
    return(x)
  }
  if(isTRUE(is.character(x))) {
    ants <- load_ants()
    return(ants$image_read(x))
  }
  if(inherits(x, "threeBrain.nii")) {
    x <- x$header
  }
  if(inherits(x, c("niftiImage", "oro.nifti"))) {
    # RNifti or oro.nifti
    tfile <- tempfile(fileext = ".nii.gz", pattern = "as_ANTsImage_")
    RNifti::writeNifti(x, file = tfile)
    on.exit({ unlink(tfile) })
    ants <- load_ants()
    return(ants$image_read(tfile))
  }
  stop("as_ANTsImage: unsupported input format: [", paste0(class(x), collapse = ","), "]. Only characters, ANTsImage, threeBrain.nii, niftiImage, or oro.nifti allowed.")
}


