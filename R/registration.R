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


#' @title Apply a transform list to map an image from one domain to another
#' @description
#' See \code{ants$apply_transforms} for more details.
#' @param fixed fixed image defining domain into which the moving image is transformed
#' @param moving moving image to be mapped to fixed space
#' @param transformlist list of strings (path to transforms) generated by
#' \code{\link{ants_registration}} where each transform is a file name
#' @param interpolator how to interpolate the image; see 'Usage'
#' @param imagetype integer: 0 (scalar), 1 (vector), 2 (tensor), 3 (time-series),
#' used when the fixed and moving images have different mode (dimensions)
#' @param whichtoinvert either \code{NULL}, \code{None} ('Python'), or a vector
#' of logical with same length as \code{transformlist}; print
#' \code{ants$apply_transforms} to see detailed descriptions
#' @param compose optional character pointing to a valid file location
#' @param defaultvalue numerical value for mappings outside the image domain
#' @param verbose whether to verbose application of transform
#' @param ... must be named arguments passing to further methods
#' @returns Transformed image. The image will share the same space as \code{fixed}.
#' @seealso \code{print(ants$apply_transforms)}
#' @examples
#'
#' if(interactive() && ants_available()) {
#'   ants <- load_ants()
#'   fixed <- as_ANTsImage( ants$get_ants_data('r16') )
#'   moving <- as_ANTsImage( ants$get_ants_data('r64') )
#'   fixed <- ants_resample_image(fixed, c(64, 64), TRUE, "linear")
#'   moving <- ants_resample_image(moving, c(64,64), TRUE, "linear")
#'
#'   mytx <- ants_registration(fixed = fixed,
#'                             moving = moving,
#'                             type_of_transform = 'SyN')
#'   mywarpedimage <- ants_apply_transforms(
#'     fixed = fixed,
#'     moving = moving,
#'     transformlist = mytx$fwdtransforms
#'   )
#'
#'   par(mfrow = c(1,3), mar = c(0,0,3,0))
#'   pal <- gray.colors(256)
#'   image(fixed[], asp = 1, axes = FALSE, col = pal,
#'         ylim = c(1, 0), main = "Reference")
#'   image(moving[], asp = 1, axes = FALSE, col = pal,
#'         ylim = c(1, 0), main = "Moving")
#'   image(mywarpedimage[], asp = 1, axes = FALSE, col = pal,
#'         ylim = c(1, 0), main = "Moving reg+resamp into Reference")
#' }
#'
#' @export
ants_apply_transforms <- function(
    fixed, moving, transformlist,
    interpolator = c("linear", "nearestNeighbor", "gaussian", "genericLabel", "bSpline",
                     "cosineWindowedSinc", "welchWindowedSinc", "hammingWindowedSinc",
                     "lanczosWindowedSinc"),
    imagetype = 0L, whichtoinvert = NULL, compose = NULL,
    defaultvalue = 0, verbose = FALSE, ...) {
  ants <- load_ants()

  verbose <- convert_if_not_python(verbose, isTRUE(as.logical(verbose)))
  defaultvalue <- convert_if_not_python(defaultvalue, {
    defaultvalue <- as.numeric(defaultvalue)
    if( is.na(defaultvalue) || defaultvalue < 0 ) { defaultvalue <- 0 }
    defaultvalue
  })


  compose <- convert_if_not_python(compose, {
    if(length(compose)) {
      compose <- normalizePath(compose, mustWork = FALSE)
    } else {
      compose <- NULL
    }
    compose
  })

  imagetype <- convert_if_not_python(imagetype, {
    imagetype <- as.integer(imagetype)
    if( length(imagetype) != 1 || is.na(imagetype) || imagetype < 0 || imagetype > 3 ) {
      stop("ants_apply_transforms: invalid `imagetype`: choose 0/1/2/3 mapping to scalar/vector/tensor/time-series")
    }
    imagetype
  })


  interpolator <- convert_if_not_python(interpolator, { match.arg(interpolator) })

  fixed <- as_ANTsImage(fixed, strict = TRUE)
  moving <- as_ANTsImage(moving, strict = TRUE)

  transformlist <- convert_if_not_python(transformlist, {
    py_list(transformlist)
  })
  n_transforms <- length(transformlist)


  if(!length(whichtoinvert)) {
    whichtoinvert <- py_none()
  } else if(length(whichtoinvert) != n_transforms){
    stop("ants_apply_transforms: `whichtoinvert` must be either NULL or have the same length as `transformlist`")
  } else {
    whichtoinvert <- py_list(as.logical(whichtoinvert))
  }

  re <- ants$apply_transforms(
    fixed = fixed,
    moving = moving,
    transformlist = transformlist,
    interpolator = interpolator,
    imagetype = imagetype,
    whichtoinvert = whichtoinvert,
    compose = compose,
    defaultvalue = defaultvalue,
    verbose = verbose,
    ...
  )
  return(re)
}



#' @title Apply a transform list to map points from one domain to another
#' @description
#' See \code{ants$apply_transforms_to_points} for more details. Please note
#' point mapping goes the opposite direction of image mapping (see
#' \code{\link{ants_apply_transforms}}), for both reasons of convention and
#' engineering.
#' @param dim dimensions of the transformation
#' @param points data frame containing columns \code{'x'}, \code{'y'},
#' \code{'z'}, \code{'t'} (depending on \code{dim})
#' @param transformlist list of strings (path to transforms) generated by
#' \code{\link{ants_registration}} where each transform is a file name
#' @param whichtoinvert either \code{NULL}, \code{None} ('Python'), or a vector
#' of logical with same length as \code{transformlist}; print
#' \code{ants$apply_transforms_to_points} to see detailed descriptions
#' @param verbose whether to verbose application of transform
#' @param ... ignored
#' @returns Transformed points in data frame (R object)
#' @seealso \code{print(ants$apply_transforms_to_points)}
#' @examples
#'
#' if(interactive() && ants_available()) {
#'   ants <- load_ants()
#'   fixed <- as_ANTsImage( ants$get_ants_data('r16') )
#'   moving <- as_ANTsImage( ants$get_ants_data('r27') )
#'
#'   reg <- ants_registration(
#'     fixed = fixed, moving = moving,
#'     type_of_transform = "antsRegistrationSyNRepro[a]")
#'
#'   pts <- data.frame(
#'     x = c(128, 127),
#'     y = c(101, 111)
#'   )
#'
#'   ptsw = ants_apply_transforms_to_points(2, pts, reg$fwdtransforms)
#'   ptsw
#' }
#'
#' @export
ants_apply_transforms_to_points <- function(
    dim, points, transformlist, whichtoinvert=NULL, verbose=FALSE, ...) {

  ants <- load_ants()

  if(inherits(dim, "python.builtin.object")) {
    dim <- py_to_r(dim)
  }
  dim <- as.integer(dim)
  verbose <- convert_if_not_python(verbose, isTRUE(as.logical(verbose)))

  if(inherits(points, "python.builtin.object")) {
    points <- py_to_r(points)
  }
  points <- as.data.frame(points)
  required_columns <- c("x", "y", "z", "t")
  required_columns <- required_columns[seq_len(dim)]
  if(!all(required_columns %in% names(points))) {
    stop("ants_apply_transforms_to_points: `points` must be a data.frame containing names: ['", paste(required_columns, collapse = "', '"), "'].")
  }
  npts <- nrow(points)
  if( npts == 0 ) { return(points) }
  if( npts == 1 ) {
    points <- rbind(points, points)
  }

  transformlist <- convert_if_not_python(transformlist, {
    py_list(transformlist)
  })
  n_transforms <- length(transformlist)


  if(!length(whichtoinvert)) {
    whichtoinvert <- py_none()
  } else if(length(whichtoinvert) != n_transforms){
    stop("ants_apply_transforms_to_points: `whichtoinvert` must be either NULL or have the same length as `transformlist`")
  } else {
    whichtoinvert <- py_list(as.logical(whichtoinvert))
  }

  re <- ants$apply_transforms_to_points(
    dim = dim,
    points = points,
    transformlist = transformlist,
    whichtoinvert = whichtoinvert,
    verbose = verbose
  )
  re <- py_to_r(re)
  if(npts == 1) {
    re <- re[1, , drop = FALSE]
  }
  return(re)
}
