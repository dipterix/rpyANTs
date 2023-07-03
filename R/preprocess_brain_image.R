#' @title Process brain image prior to segmentation
#' @description
#' Strip skulls, normalize intensity, align and re-sample to template. This
#' procedure is needed for many \code{antspynet} functions since the deep
#' neural networks are trained in template spaces
#' @param x \code{'ANTsImage'} or path to image to process
#' @param truncate_intensity defines the quantile threshold for truncating the
#' image intensity
#' @param brain_extraction_modality character of length 1, perform brain
#' extraction modality
#' @param template_transform_type either \code{'Rigid'} or \code{'Affine'}
#' align to template brain
#' @param template template image (not skull-stripped) or string, e.g.
#' \code{'biobank'}, \code{'croppedMni152'}
#' @param do_bias_correction whether to perform bias field correction
#' @param return_bias_field return bias field as an additional output without
#' bias correcting the image
#' @param do_denoising whether to remove noises using non-local means
#' @param intensity_matching_type either \code{'regression'} or
#' \code{'histogram'}; only is performed if \code{reference_image} is not
#' \code{NULL}.
#' @param reference_image \code{'ANTsImage'} or path to image, or \code{NULL}
#' @param intensity_normalization_type either re-scale the intensities to
#' \code{c(0, 1)} (\code{'01'}), or for zero-mean, unit variance
#' (\code{'0mean'}); if \code{NULL} normalization is not performed
#' @param verbose print progress to the screen
#' @returns Dictionary with images after process. The images are registered and
#' re-sampled into template.
#' @seealso \code{antspynet$preprocess_brain_image}
#' @examples
#'
#' library(rpyANTs)
#' if(interactive() && ants_available("antspynet")) {
#'   image_path <- ants$get_ants_data('r30')
#'   preprocessed <- antspynet_preprocess_brain_image(
#'     image_path, verbose = FALSE
#'   )
#'
#'   # Compare
#'   orig_img <- as_ANTsImage(image_path)
#'   new_img <- preprocessed$preprocessed_image
#'   pal <- grDevices::gray.colors(256, start = 0, end = 1)
#'
#'   par(mfrow = c(1, 2), mar = c(0.1, 0.1, 0.1, 0.1),
#'       bg = "black", fg = "white")
#'   image(orig_img[], asp = 1, axes = FALSE,
#'         col = pal, ylim = c(1, 0))
#'   image(new_img[], asp = 1, axes = FALSE,
#'         col = pal, ylim = c(1, 0))
#'
#' }
#'
#'
#' @export
antspynet_preprocess_brain_image <- function(
    x, truncate_intensity = c(0.01, 0.99),
    brain_extraction_modality = c(
      "none", "t1", "t1v0", "t1nobrainer", "t1combined", "flair", "t2",
      "bold", "fa", "t1infant", "t2infant"),
    template_transform_type = c("None", "Affine", "Rigid"),
    template = c("biobank", "croppedMni152"),
    do_bias_correction = TRUE,
    return_bias_field = FALSE,
    do_denoising = TRUE,
    intensity_matching_type = c("regression", "histogram"),
    reference_image = NULL,
    intensity_normalization_type = NULL,
    verbose = TRUE
) {

  brain_extraction_modality <- convert_if_not_python(
    brain_extraction_modality,
    {
      brain_extraction_modality <- match.arg(brain_extraction_modality)
      if(brain_extraction_modality == "none") {
        brain_extraction_modality <- NULL
      }
      brain_extraction_modality
    }
  )

  template_transform_type <- convert_if_not_python(
    template_transform_type,
    {
      if(length(template_transform_type) > 1) {
        template_transform_type <- template_transform_type[[1]]
      }
      if(template_transform_type == "None") {
        template_transform_type <- NULL
      }
      template_transform_type
    }
  )

  intensity_matching_type <- convert_if_not_python(
    intensity_matching_type,
    {
      match.arg(intensity_matching_type)
    }
  )

  if(missing(template) || (is.character(template) && length(template) > 1 )) {
    template <- match.arg(template)
  }

  do_bias_correction <- convert_if_not_python(
    do_bias_correction, isTRUE(as.logical(do_bias_correction)))

  verbose <- convert_if_not_python(
    verbose, isTRUE(as.logical(verbose)))

  return_bias_field <- convert_if_not_python(
    return_bias_field, isTRUE(as.logical(return_bias_field)))

  do_denoising <- convert_if_not_python(
    do_denoising, isTRUE(as.logical(do_denoising)))

  truncate_intensity <- convert_if_not_python(
    truncate_intensity, do.call(tuple, unname(as.list(truncate_intensity)))
  )

  intensity_normalization_type <- convert_if_not_python(
    intensity_normalization_type, {
      if(length(intensity_normalization_type)) {
        as.character(intensity_normalization_type)
      } else {
        NULL
      }
    })

  antspynet <- load_antspynet()
  image_ants <- as_ANTsImage(x, strict = TRUE)
  reference_image <- as_ANTsImage(reference_image, strict = FALSE)

  re <- antspynet$preprocess_brain_image(
    image = image_ants,
    truncate_intensity = truncate_intensity,
    brain_extraction_modality = brain_extraction_modality,
    template_transform_type = template_transform_type,
    template = template,
    do_bias_correction = do_bias_correction,
    return_bias_field = return_bias_field,
    do_denoising = do_denoising,
    intensity_matching_type = intensity_matching_type,
    reference_image = reference_image,
    intensity_normalization_type = intensity_normalization_type,
    verbose = verbose
  )
  return(re)
}
