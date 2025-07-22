#' @title Extract brain and strip skull
#' @name brain_extraction
#' @description
#' Function \code{brain_mask} and \code{brain_extraction} use the \code{ants}
#' base package.
#' Function \code{antspynet_brain_extraction} uses \code{antspynet} to extract
#' the brain using deep neural network. This requires additional configuration.
#' Print \code{antspynet$brain_extraction} to see the original documentation.
#' @param x input image or image path
#' @param work_path working directory; default is temporary path
#' @param auto_clean whether to automatically clean the working path if
#' the path is temporary; default is true
#' @param skull_alpha used by \code{brain_extraction}, the opacity of the
#' skull in the final image; default is 0 (completely strip the skull);
#' set to values between 0 to 1 to add skulls (in such case, the background
#' noises will be set to 0)
#' @param threshold_quantile used only when \code{skull_alpha} is positive to
#' remove the noises
#' @param modality modality type, used by \code{antspynet_brain_extraction} only
#' @param verbose whether to print out process to the screen
#' @param ... see \code{work_path} and \code{auto_clean}
#' @returns \describe{
#' \item{\code{brain_mask}}{Brain mask image}
#' \item{\code{brain_extraction}}{extracted brain}
#' \item{\code{antspynet_brain_extraction}}{Brain mask image}
#' }
#' @export
antspynet_brain_extraction <- function(
    x, modality = c("t1", "t1nobrainer", "t1combined", "flair", "t2", "t2star", "bold", "fa", "t1t2infant", "t1infant", "t2infant"),
    verbose=FALSE) {

  antspynet <- load_antspynet()
  image <- as_ANTsImage(x)
  modality <- convert_if_not_python(modality, { match.arg(modality) })
  verbose <- convert_if_not_python(verbose, { isTRUE(as.logical(verbose)) })
  re <- antspynet$brain_extraction(image = image, modality = modality, verbose = verbose)
  return(re)
}

#' @rdname brain_extraction
#' @export
brain_mask <- function(
    x, work_path = tempfile(pattern = "rpyant_brain_extraction_"),
    verbose = TRUE, auto_clean = TRUE) {

  # x <- "/Users/dipterix/rave_data/raw_dir/PAV058/rave-imaging/inputs/MRI/MRI_RAW.nii"
  # work_path <- tempfile(pattern = "rpyant_brain_extraction_")
  # use_temppath <- TRUE
  # verbose <- TRUE
  # auto_clean <- FALSE

  use_temppath <- missing(work_path)

  # MNI152a should do the work
  template_path <- ensure_template("mni_icbm152_nlin_asym_09a")
  ants <- load_ants()
  rpyants <- load_rpyants()

  if(use_temppath && auto_clean) {
    dir_create2(work_path, showWarnings = FALSE, recursive = TRUE)
    on.exit({
      unlink(work_path, recursive = TRUE)
    })
  }
  work_path <- normalize_path(work_path, must_work = TRUE)
  template_path <- normalize_path(template_path, must_work = TRUE)

  if( verbose ) {
    verbose <- TRUE
  } else {
    verbose <- FALSE
  }

  template_t1w <- file_path(template_path, "T1.nii.gz")
  template_mask <- file_path(template_path, "T1_brainmask.nii.gz")
  image <- as_ANTsImage(x, strict = TRUE)


  mask <- rpyants$common$brainmask(
    # img_path
    image,

    # template_path,
    template_t1w,

    template_mask,

    # working_path,
    work_path,

    # verbose
    verbose
  )

  return(mask)
}

#' @rdname brain_extraction
#' @export
brain_extraction <- function(x, skull_alpha = 0,
                             threshold_quantile = 0.85, ..., verbose = TRUE) {
  image <- as_ANTsImage(x, strict = TRUE)
  mask <- brain_mask(x = image, ..., verbose = verbose)

  brain <- image * mask

  if( skull_alpha > 0 ) {
    if( skull_alpha > 1 ) {
      skull_alpha <- 1
    }
    pixel_values <- image[mask == 0]
    pixel_values <- pixel_values[pixel_values > 0]
    threshold <- 0
    if(length(pixel_values)) {
      threshold_quantile <- as.numeric(threshold_quantile)
      if(threshold_quantile < 0) {
        threshold_quantile <- 0
      } else if (threshold_quantile > 1) {
        threshold_quantile <- 1
      }
      threshold <- stats::quantile(pixel_values, threshold_quantile)

      skull_mask <- (image > threshold)
      skull_mask[mask] <- 0
      skull <- image * skull_mask * skull_alpha
      brain <- brain + skull
    }

  }

  return(brain)
}
