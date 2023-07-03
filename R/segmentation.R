#' @name antspynet_segmentation
#' @title Imaging segmentation using \code{antspynet}
#' @description
#' Supports \code{Desikan-Killiany-Tourville} labeling and deep \code{'Atropos'}.
#' @param x \code{'NIfTI'} image or path to the image that is to be segmented
#' @param do_preprocessing whether \code{x} is in native space and needs the be
#' registered to template brain before performing segmentation; default is true
#' since the model is trained with template brain. If you want to manually
#' process the image, see \code{\link{antspynet_preprocess_brain_image}}
#' @param return_probability_images whether to return probability images
#' @param do_lobar_parcellation whether to perform lobar 'parcellation'
#' @param use_spatial_priors whether to use \code{'MNI'} partial tissue priors
#' @param aseg_only whether to just return the segmented image
#' @param verbose whether to print out the messages
#' @return One or a list of \code{'ANTsImage'} image instances. Please print
#' out \code{antspynet$desikan_killiany_tourville_labeling} or
#' \code{antspynet$deep_atropos} to see the details.
#' @seealso \code{antspynet$desikan_killiany_tourville_labeling},
#' \code{antspynet$deep_atropos}
#' @examples
#'
#'
#' # Print Python documents
#' if(interactive() && ants_available("antspynet")) {
#'   antspynet <- load_antspynet()
#'
#'   print(antspynet$deep_atropos)
#'
#'   print(antspynet$desikan_killiany_tourville_labeling)
#' }
#'
#'
NULL

#' @rdname antspynet_segmentation
#' @export
antspynet_desikan_killiany_tourville_labeling <- function(
    x, do_preprocessing = TRUE, return_probability_images = FALSE,
    do_lobar_parcellation = FALSE, verbose = TRUE) {
  # DIPSAUS DEBUG START
  # x <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/raw/Suthana_Demo/T1.nii.gz"

  antspynet <- load_antspynet()

  do_preprocessing <- convert_if_not_python(do_preprocessing, {
    isTRUE(as.logical(do_preprocessing))
  })

  return_probability_images <- convert_if_not_python(return_probability_images, {
    isTRUE(as.logical(return_probability_images))
  })

  do_lobar_parcellation <- convert_if_not_python(do_lobar_parcellation, {
    isTRUE(as.logical(do_lobar_parcellation))
  })

  verbose <- convert_if_not_python(verbose, {
    isTRUE(as.logical(verbose))
  })

  image_ants <- as_ANTsImage(x, strict = TRUE)

  # Performing AntsPyNet Segmentation
  antspynet_segmentation <- antspynet$desikan_killiany_tourville_labeling(
    t1 = image_ants, do_preprocessing = do_preprocessing,
    return_probability_images = return_probability_images,
    do_lobar_parcellation = do_lobar_parcellation, verbose = verbose)

  return(antspynet_segmentation)
}

#' @rdname antspynet_segmentation
#' @export
antspynet_deep_atropos <- function(
    x, do_preprocessing = TRUE, use_spatial_priors = TRUE,
    aseg_only = TRUE, verbose = TRUE) {
  # DIPSAUS DEBUG START
  # x <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/raw/Suthana_Demo/T1.nii.gz"

  antspynet <- load_antspynet()
  image_ants <- as_ANTsImage(x, strict = TRUE)

  do_preprocessing <- convert_if_not_python(do_preprocessing, {
    isTRUE(as.logical(do_preprocessing))
  })

  use_spatial_priors <- convert_if_not_python(use_spatial_priors, {
    use_spatial_priors <- as.integer(use_spatial_priors)
    if( is.na(use_spatial_priors) || use_spatial_priors < 0 || use_spatial_priors > 1) {
      use_spatial_priors <- 1L
    }
    use_spatial_priors
  })

  verbose <- convert_if_not_python(verbose, {
    isTRUE(as.logical(verbose))
  })



  # Performing AntsPyNet Segmentation
  antspynet_segmentation <- antspynet$deep_atropos(image_ants, do_preprocessing = do_preprocessing, use_spatial_priors = use_spatial_priors, verbose = verbose)

  if( aseg_only ) {
    antspynet_segmentation <- antspynet_segmentation$segmentation_image
  }

  return(antspynet_segmentation)
}

