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
#'
py_slice <- function(...) {
  reticulate::import_builtins()$slice(...)
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


