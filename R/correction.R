#' Truncate and correct 'MRI' intensity
#' @description
#' Uses \code{ants.abp_n4} to truncate and correct intensity
#' @param image 'MRI' image to be corrected, will be passed to
#' \code{\link{as_ANTsImage}}
#' @param mask binary mask image
#' @param intensity_truncation numerical length of two, quantile probabilities
#' to truncate.
#' @returns An \code{'ANTsImage'} instance
#' @examples
#'
#' if(interactive() && ants_available()) {
#'   ants <- load_ants()
#'   scale <- (0.1 + outer(
#'     seq(0, 1, length.out = 256)^6,
#'     seq(0, 1, length.out = 256)^2,
#'     FUN = "+"
#'   )) / 6
#'   img = ants$image_read(ants$get_ants_data('r16')) * scale
#'
#'   corrected <- correct_intensity(img)
#'
#'   pal <- gray.colors(255, start = 0)
#'   par(mfrow = c(1, 2), mar = c(0.1, 0.1, 2.1, 0.1),
#'       bg = "black", fg = "white")
#'   image(img[], asp = 1, axes = FALSE,
#'         col = pal, ylim = c(1, 0),
#'         main = "Original", col.main = "white")
#'   image(corrected[], asp = 1, axes = FALSE,
#'         col = pal, ylim = c(1, 0),
#'         main = "Corrected", col.main = "white")
#'
#' }
#'
#' @export
correct_intensity <- function(
  image, mask = NULL, intensity_truncation = c(0.025, 0.975)
) {
  # DIPSAUS DEBUG START
  # image <- "/Users/dipterix/rave_data/raw_dir/testtest2/rave-imaging/inputs/anat/sub-testtest2_ses-preop_acq-ax_desc-preproc_T1w.nii.gz"
  # intensity_truncation = c(0.025, 0.975)
  # as_uint8 = TRUE
  # mask <- "/Users/dipterix/rave_data/raw_dir/testtest2/rave-imaging/ants/mri/brainmask.nii.gz"
  assert_quantile2 <- function(qt) {
    stopifnot(
      length(qt) == 2 && is.numeric(qt) &&
        isTRUE(qt[[1]] < qt[[2]]) &&
        qt[[1]] < 1 && qt[[1]] >= 0 &&
        qt[[2]] <= 1 && qt[[2]] > 0
    )
  }

  intensity_truncation <- as.double(intensity_truncation)
  assert_quantile2(intensity_truncation)
  image <- as_ANTsImage(image, strict = TRUE)
  mask <- as_ANTsImage(mask, strict = FALSE)

  # N4 bias correction
  res <- image$abp_n4(
    mask = mask,
    intensity_truncation = list(intensity_truncation[[1]],
                                intensity_truncation[[2]], 256L)
  )
  res
}
