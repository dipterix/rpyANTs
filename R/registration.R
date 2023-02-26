#' @title Register two images using 'ANTs'
#' @param fixed fixed image to which we register the moving image, can be
#' character path to 'NIfTI' image, or \code{'ANTsImage'} instance,
#' \code{'oro.nifti'} object, \code{'niftiImage'} from
#' package \code{'RNifti'}, or \code{'threeBrain.nii'} from package
#' \code{'threeBrain'}; see also \code{\link{as_ANTsImage}}
#' @param moving moving image to be mapped to fixed space; see also \code{\link{as_ANTsImage}}
#' @param type_of_transform a linear or non-linear registration type;
#' print \code{ants$registration} to see details
#' @param initial_transform optional list of strings; transforms to apply prior
#' to registration
#' @param outprefix output file to save results
#' @param mask image mask; see also \code{\link{as_ANTsImage}}
#' @param grad_step,flow_sigma,total_sigma optimization parameters
#' @param aff_metric the metric for the 'affine' transformation, choices are
#' \code{'GC'}, \code{'mattes'}, \code{'meansquares'}
#' @param aff_sampling,aff_random_sampling_rate,aff_iterations,aff_shrink_factors,aff_smoothing_sigmas controls 'affine' transform
#' @param syn_metric the metric for the \code{'SyN'} transformation, choices
#' are \code{'GC'}, \code{'mattes'}, \code{'meansquares'}, \code{'demons'}
#' @param syn_sampling,reg_iterations controls the \code{'SyN'} transform
#' @param write_composite_transform whether the composite transform (and its
#' inverse, if it exists) should be written to an 'HDF5' composite file;
#' default is false
#' @param verbose verbose the progress
#' @param smoothing_in_mm logical, currently only impacts low dimensional
#' registration
#' @param ... others passed to \code{ants$registration}
#'
#' @returns A 'Python' dictionary of aligned images and transform files.
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   ants <- load_ants()
#'
#'   # check the python documentation here for detailed explanation
#'   print(ants$registration)
#'
#'   # example to register
#'   fi <- ants$image_read(ants$get_ants_data('r16'))
#'   mo <- ants$image_read(ants$get_ants_data('r64'))
#'
#'   # resample to speed up this example
#'   fi <- ants$resample_image(fi, list(60L,60L), TRUE, 0L)
#'   mo <- ants$resample_image(mo, list(60L,60L), TRUE, 0L)
#'
#'   # SDR transform
#'   transform <- ants_registration(
#'     fixed=fi, moving=mo, type_of_transform = 'SyN' )
#'
#'   ants$plot(fi, overlay = transform$warpedmovout, overlay_alpha = 0.3)
#'
#'
#' }
#'
#'
#' @export
ants_registration <- function(
    fixed, moving, type_of_transform='SyN', initial_transform = NULL,
    outprefix=tempfile(), mask=NULL, grad_step=0.2, flow_sigma=3, total_sigma=0,
    aff_metric=c('mattes', 'GC', 'meansquares'),
    aff_sampling=32, aff_random_sampling_rate=0.2,
    syn_metric=c('mattes', 'CC', 'meansquares', 'demons'),
    syn_sampling=32, reg_iterations=c(40, 20, 0),
    aff_iterations=c(2100, 1200, 1200, 10),
    aff_shrink_factors=c(6, 4, 2, 1),
    aff_smoothing_sigmas=c(3, 2, 1, 0),
    write_composite_transform=FALSE, verbose=FALSE,
    smoothing_in_mm=FALSE, ...) {

  # DIPSAUS DEBUG START
  # moving <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/raw/PAV006/rave-imaging/derivative/CT_RAW.nii.gz"
  # fixed <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/raw/PAV006/rave-imaging/derivative/MRI_RAW.nii"
  # type_of_transform <- "Rigid"

  aff_metric <- convert_if_not_python(aff_metric, { match.arg(aff_metric) })
  syn_metric <- convert_if_not_python(syn_metric, { match.arg(syn_metric) })

  grad_step <- convert_if_not_python(grad_step, as.double(grad_step))
  flow_sigma <- convert_if_not_python(flow_sigma, as.double(flow_sigma))
  total_sigma <- convert_if_not_python(total_sigma, as.double(total_sigma))
  aff_sampling <- convert_if_not_python(aff_sampling, as.integer(aff_sampling))
  aff_random_sampling_rate <- convert_if_not_python(aff_random_sampling_rate, as.double(aff_random_sampling_rate))
  syn_sampling <- convert_if_not_python(syn_sampling, as.integer(syn_sampling))

  reg_iterations <- convert_if_not_python(reg_iterations, tuple(as.list(as.integer(reg_iterations))))
  aff_iterations <- convert_if_not_python(aff_iterations, tuple(as.list(as.integer(aff_iterations))))
  aff_shrink_factors <- convert_if_not_python(aff_shrink_factors, tuple(as.list(as.integer(aff_shrink_factors))))
  aff_smoothing_sigmas <- convert_if_not_python(aff_smoothing_sigmas, tuple(as.list(as.integer(aff_smoothing_sigmas))))

  write_composite_transform <- convert_if_not_python(write_composite_transform, as.logical(write_composite_transform))
  verbose <- convert_if_not_python(verbose, as.logical(verbose))
  smoothing_in_mm <- convert_if_not_python(smoothing_in_mm, as.logical(smoothing_in_mm))


  type_of_transform <- convert_if_not_python(type_of_transform, as.character(type_of_transform))

  fixed_img <- as_ANTsImage(fixed, strict = TRUE)
  moving_img <- as_ANTsImage(moving, strict = TRUE)
  mask <- as_ANTsImage(mask, strict = FALSE)

  if(length(initial_transform)) {
    initial_transform <- as_ANTsTransform(initial_transform, fixed_img$dimension)
  } else {
    initial_transform <- NULL
  }

  outprefix <- convert_if_not_python(outprefix, {
    normalizePath(outprefix, mustWork = FALSE, winslash = "/")
  })


  ants <- load_ants()
  tfiles1 <- snapshot_tempfiles()
  py_results <- ants$registration(
    fixed = fixed_img, moving = moving_img,
    type_of_transform = type_of_transform,
    initial_transform = initial_transform,
    outprefix = outprefix, mask = mask,
    grad_step=grad_step, flow_sigma=flow_sigma, total_sigma=total_sigma,
    aff_metric=aff_metric, aff_sampling=aff_sampling,
    aff_random_sampling_rate=aff_random_sampling_rate,
    syn_metric=syn_metric, syn_sampling=syn_sampling,
    reg_iterations=reg_iterations,
    aff_iterations=aff_iterations,
    aff_shrink_factors=aff_shrink_factors,
    aff_smoothing_sigmas=aff_smoothing_sigmas,
    write_composite_transform=write_composite_transform,
    verbose=verbose,
    smoothing_in_mm=smoothing_in_mm)
  tfiles2 <- snapshot_tempfiles()
  remove_tmpfiles(setdiff(tfiles2, tfiles1))

  py_results

}
