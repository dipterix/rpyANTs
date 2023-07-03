#' @title Resample image
#' @description
#' See \code{ants$resample_image} for more details
#' @param x input image
#' @param resample_params either relative number or absolute integers
#' @param use_voxels whether the \code{resample_params} should be treated as
#' new dimension \code{use_voxels=TRUE}, or the new dimension should be
#' calculated based on current dimension and \code{resample_params} combined
#' (\code{use_voxels=FALSE} then \code{resample_params} will be treated as
#' relative number); default is \code{FALSE}
#' @param interp_type interpolation type; either integer or character; see
#' 'Usage' for available options
#' @returns Resampled image
#' @examples
#'
#'
#' if(interactive() && ants_available()) {
#'   ants <- load_ants()
#'   fi <- as_ANTsImage(ants$get_ants_data("r16"))
#'
#'   # linear (interp_type = 0 or "linear)
#'   filin <- ants_resample_image(fi, c(50, 60), TRUE, "linear")
#'
#'   # nearest neighbor (interp_type = 1 or "nn)
#'   finn <- ants_resample_image(fi, c(50, 60), TRUE, "nn")
#'
#'   par(mfrow = c(1, 3), mar = c(0, 0, 0, 0))
#'   pal <- gray.colors(256, start = 0)
#'
#'   image(fi[], asp = 1, axes = FALSE,
#'         ylim = c(1,0), col = pal)
#'   image(filin[], asp = 1, axes = FALSE,
#'         ylim = c(1,0), col = pal)
#'   image(finn[], asp = 1, axes = FALSE,
#'         ylim = c(1,0), col = pal)
#' }
#'
#'
#' @export
ants_resample_image <- function(x, resample_params, use_voxels = FALSE,
                                interp_type = c("linear", "nn", "guassian", "sinc", "bspline")) {
  # DIPSAUS DEBUG START
  # x <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/raw/Suthana_Demo/T1.nii.gz"
  # type <- "bspline"

  interp_type2 <- as.character(interp_type)
  if(!length(interp_type2)) {
    interp_type2 <- "linear"
  } else if(length(interp_type2) > 1 || !grepl("^[1-5]", interp_type2) ) {
    interp_type2 <- match.arg(arg = interp_type2,
                              choices = c("linear", "nn", "guassian", "sinc", "bspline"),
                              several.ok = FALSE)
  } else {
    interp_type2 <- as.integer(interp_type2)
  }
  if(is.na(interp_type2) || !is.integer(interp_type2)) {
    interp_type2 <- which(c("linear", "nn", "guassian", "sinc", "bspline") == interp_type2)
  }
  if(is.na(interp_type2)) {
    stop("ants_resample_image: Invalid `interp_type`: ", as.character(interp_type))
  }

  resample_params <- convert_if_not_python(resample_params, {
    if(!is.numeric(resample_params) || any(resample_params < 0)) {
      stop("ants_resample_image: Invalid `resample_params`: all dimensions must be positive")
    }
    as.numeric(resample_params)
  })

  use_voxels <- convert_if_not_python(use_voxels, {
    isTRUE(as.logical(use_voxels))
  })

  image_ants <- as_ANTsImage(x, strict = TRUE)

  re <- image_ants$resample_image(
    resample_params = resample_params,
    use_voxels = use_voxels,
    interp_type = interp_type2
  )

  return(re)
}
