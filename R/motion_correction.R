#' Motion correction
#' @description
#' Print \code{ants$motion_correction} to see the original document
#'
#' @param x input image, usually 'fMRI' series
#' @param fixed fixed image to register all timepoints to
#' @param type_of_transform see \code{\link{ants_registration}}
#' @param mask mask for image
#' @param fdOffset offset value to use in frame-wise displacement calculation
#' @param outprefix save path
#' @param verbose whether to verbose the messages
#' @param ... passed to registration methods
#' @returns Motion-corrected image
#' @examples
#'
#' if(interactive() && ants_available()) {
#'   fi <- as_ANTsImage(ants$get_ants_data('ch2'))
#'   mytx <- ants_motion_correction( fi )
#'
#'   par(mfrow = c(1, 2), mar = c(1,1,1,1))
#'   image(fi[,,91], asp = 1, axes = FALSE)
#'   image(mytx$motion_corrected[,,91], asp = 1, axes = FALSE)
#' }
#'
#' @export
ants_motion_correction <- function(
    x,
    fixed = NULL,
    type_of_transform = 'BOLDRigid',
    mask = NULL,
    fdOffset = 50,
    outprefix = '',
    verbose = FALSE,
    ...
) {

  ants <- load_ants()

  x <- as_ANTsImage(x, strict = TRUE)
  fixed <- as_ANTsImage(fixed, strict = FALSE)
  mask <- as_ANTsImage(mask, strict = FALSE)
  fdOffset <- convert_if_not_python(fdOffset, as.numeric(fdOffset))
  outprefix <- convert_if_not_python(outprefix, {
    outprefix_new <- ""
    if(length(outprefix) >= 1) {
      outprefix <- outprefix[[1]]
      if(!is.na(outprefix) && trimws(outprefix) != "") {
        outprefix_new <- normalizePath(outprefix, winslash = "/", mustWork = FALSE)
      }
    }
    outprefix_new
  })

  verbose <- convert_if_not_python(verbose, {
    isTRUE(as.logical(verbose))
  })

  type_of_transform <- convert_if_not_python(type_of_transform, {
    type_of_transform <- as.character(type_of_transform)
    if(length(type_of_transform) >= 1) {
      type_of_transform <- type_of_transform[[1]]
    } else if(length(type_of_transform) == 0) {
      type_of_transform <- 'BOLDRigid'
    }
    type_of_transform
  })


  ants$motion_correction(
    image = x,
    fixed = fixed,
    mask = mask,
    type_of_transform = type_of_transform,
    fdOffset = fdOffset,
    outprefix = outprefix,
    verbose = verbose,
    ...
  )
}
