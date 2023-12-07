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
#' @details
#' Function family \code{ants_registration*} align images (specified by
#' \code{moving}) to \code{fixed}. Here are descriptions of the variations:
#' \describe{
#' \item{\code{ants_registration}}{Simple wrapper function for 'Python'
#' implementation \code{ants.registration}, providing various of registration
#' options}
#' \item{\code{ants_registration_halpern1}}{Rigid-body registration designed
#' for 'Casey-Halpern' lab, mainly used for aligning 'MRI' to 'CT' (or the other
#' way around)}
#' }
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
    if(length(initial_transform) == 1 && is.character(initial_transform) && file.exists(initial_transform)) {
      initial_transform <- as_ANTsTransform(initial_transform, fixed_img$dimension)
    }
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

#' @name halpern_preprocess
#' @title 'ANTs' functions for 'Halpern' lab
#' @param fixed fixed image as template
#' @param moving moving image that is to be registered into \code{fixed}
#' @param outprefix output prefix, needs to be absolute path prefix
#' @param verbose whether to verbose the progress; default is true
#' @param fixed_is_ct whether \code{fixed} is 'CT'
#' @param mask mask file for template (skull-stripped)
#' @param roi_folder template 'ROI' or atlas folder in which the image atlases
#' or masks will be transformed into subject's native brain
#' @returns A list of result configurations
#' @export
halpern_register_ct_mri <- function(fixed, moving, outprefix, fixed_is_ct = TRUE, verbose = TRUE) {
  # DIPSAUS DEBUG START
  # moving <- "~/rave_data/raw_dir/testtest/rave-imaging/coregistration/MRI_RAW.nii.gz"
  # fixed <- "~/rave_data/raw_dir/testtest/rave-imaging/coregistration/CT_RAW.nii.gz"
  # outprefix <- "~/rave_data/raw_dir/testtest/rave-imaging/coregistration/t1w_to_postopct_"
  # verbose = TRUE
  # fixed_is_ct <- TRUE

  rpyants_py <- load_rpyants()

  outprefix <- normalize_path(outprefix, must_work = FALSE)
  dir.create(dirname(outprefix), showWarnings = FALSE, recursive = TRUE)
  verbose <- convert_if_not_python(verbose, as.logical(verbose))

  fixed_img <- as_ANTsImage(fixed, strict = TRUE)
  moving_img <- as_ANTsImage(moving, strict = TRUE)

  fixed <- sprintf("%sorig_fixed.nii.gz", outprefix)
  moving <- sprintf("%sorig_moving.nii.gz", outprefix)

  fixed_img$to_file(filename = fixed)
  moving_img$to_file(filename = moving)

  py_results <- rpyants_py$registration$halpern_coregister_ct_mri(
    fixed_path = fixed,
    moving_path = moving,
    outprefix = outprefix
  )

  transform <- py_to_r(py_results$transforms)
  ct_lps_to_mri_lps <- as.matrix(as_ANTsTransform(transform))
  if( fixed_is_ct ) {
    fixed_is_ct <- TRUE
    ct_ijk_to_lps <- t(t(py_to_r(fixed_img$direction)) *
                         as.double(py_to_r(fixed_img$spacing)))
    ct_ijk_to_lps <- rbind(cbind(ct_ijk_to_lps, as.double(py_to_r(fixed_img$origin))), c(0, 0, 0, 1))
  } else {
    fixed_is_ct <- FALSE
    ct_lps_to_mri_lps <- solve(ct_lps_to_mri_lps)
    ct_ijk_to_lps <- t(t(py_to_r(moving_img$direction)) *
                         as.double(py_to_r(moving_img$spacing)))
    ct_ijk_to_lps <- rbind(cbind(ct_ijk_to_lps, as.double(py_to_r(moving_img$origin))), c(0, 0, 0, 1))
  }

  ct_ijk_to_mri_lps <- ct_lps_to_mri_lps %*% ct_ijk_to_lps
  ct_ijk_to_mri_ras <- diag(c(-1, -1, 1, 1)) %*% ct_ijk_to_mri_lps

  utils::write.table(ct_ijk_to_mri_ras, paste0(outprefix, "CT_IJK_to_MR_RAS.txt"),
                     row.names = FALSE, col.names = FALSE)

  ct_ras_to_mri_ras <- diag(c(-1, -1, 1, 1)) %*%
    ct_lps_to_mri_lps %*% diag(c(-1, -1, 1, 1))
  utils::write.table(ct_ras_to_mri_ras, paste0(outprefix, "CT_RAS_to_MR_RAS.txt"),
                     row.names = FALSE, col.names = FALSE)

  # convert to R object
  structure(
    list(
      description = "Rigid-body CT-MRI co-registration",
      prefix = outprefix,
      fixed_is_ct = fixed_is_ct,
      fixed = "orig_fixed.nii.gz",
      moving = "orig_moving.nii.gz",
      warped = "Warped.nii.gz",
      transform = list(
        raw = "0GenericAffine.mat",
        ct_vox_to_mri_ras = "CT_IJK_to_MR_RAS.txt",
        ct_ras_to_mri_ras = "CT_RAS_to_MR_RAS.txt"
      )
    ),
    class = "rpyANTs.halpern_register_ct_mri"
  )
}

#' @rdname halpern_preprocess
#' @export
halpern_register_template_mri <- function(fixed, moving, outprefix, mask = NULL, verbose = TRUE) {

  # DIPSAUS DEBUG START
  # moving <- "~/rave_data/raw_dir/testtest/rave-imaging/coregistration/MRI_RAW.nii.gz"
  # fixed <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/Tower-related/templates 2/Lead-DBS_atlases_MNI_ICBM_2009b_NLIN_ASYM/MNI_ICBM_2009b_NLIN_ASYM/t1_brain.nii.gz"
  # outprefix = "/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/outputs/t1_to_mri_"
  # verbose <- TRUE
  # type_of_transform <- "Rigid"
  # mask <- NULL
  # raveio::backup_file(dirname(outprefix), remove = TRUE)
  # mask <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/Tower-related/templates 2/Lead-DBS_atlases_MNI_ICBM_2009b_NLIN_ASYM/MNI_ICBM_2009b_NLIN_ASYM/t1_brain_mask.nii.gz"

  ants <- load_ants()

  verbose <- convert_if_not_python(verbose, as.logical(verbose))

  fixed_img <- as_ANTsImage(fixed, strict = TRUE)
  moving_img <- as_ANTsImage(moving, strict = TRUE)
  mask_img <- as_ANTsImage(mask, strict = FALSE)

  # # initial affine register
  # AffBasename = 't1_to_MNI_aff_'
  # OutputAff = str(dpReg + AffBasename)
  # ANTs_command = AntsPth + str('antsRegistrationSyN.sh -d 3 -m %fnBrainNii -f %fnMNI -o %OutputAff -n %numcores -t a -r 2 -j 1') % (fnBrainNii, fnMNI, OutputAff, numcores)
  #
  # SyNbasename = 't1_to_MNI_opt_';
  # AffBrainNii = str(dpReg + str(str(AffBasename) + 'Warped.nii.gz'))
  # OutNonl = str(dpReg + SyNbasename)
  # cmd = AntsPth + str('antsRegistrationSyN.sh -d 3 -m %s -f %s -o %s -n %s -t b -r 2 -j 1 -x %s') % (AffBrainNii, fnMNI, OutNonl, numcores, fnVDCmask)

  outprefix <- normalizePath(to_r(outprefix), mustWork = FALSE, winslash = "/")
  dir.create(dirname(outprefix), showWarnings = FALSE, recursive = TRUE)

  outprefix1 <- r_to_py(sprintf("%saffine_", outprefix))
  outprefix2 <- r_to_py(sprintf("%sdeformable_", outprefix))

  # -t a -r 2 -j 1
  affine_reg <- ants$registration(
    fixed = fixed_img,
    moving = moving_img,
    mask = mask_img,
    outprefix = outprefix1,
    type_of_transform = "antsRegistrationSyN[a]",
    aff_metric = "CC",
    syn_metric = "CC",
    write_composite_transform = FALSE,
    aff_sampling = 2L,
    syn_sampling = 2L,
    verbose = verbose
  )

  affine_reg_r <- py_to_r(affine_reg)
  # Save transformed outputs
  affine_reg_r$warpedmovout$to_file(sprintf("%saffine_warpedmovout.nii.gz", outprefix))
  fixed_img$to_file(sprintf("%sorig_fixed.nii.gz", outprefix))
  moving_img$to_file(sprintf("%sorig_moving.nii.gz", outprefix))

  # -t b -r 2 -j 1 -x mask
  syn_reg <- ants$registration(
    fixed = fixed_img,
    moving = affine_reg_r$warpedmovout,
    mask = mask_img,
    outprefix = outprefix2,
    type_of_transform = "antsRegistrationSyN[b]",
    aff_metric = "CC",
    syn_metric = "CC",
    write_composite_transform = FALSE,
    aff_sampling = 2L,
    syn_sampling = 2L,
    verbose = verbose,
  )
  syn_reg_r <- py_to_r(syn_reg)

  # Save transformed outputs
  syn_reg_r$warpedmovout$to_file(sprintf("%sdeformable_warpedmovout.nii.gz", outprefix))

  # 'antsApplyTransforms -d 3 -r %s -i \'%s\' -o %s -n NearestNeighbor -t [%s,1] -t [%s,1] -t [%s,1] -t %s' % (fnEpiBrainNii, fnROImask[i], fnOutput, fnBbrMat, fnAff_MNIaff, fnSyN_MNIaff, fnSyN_MNIwarp)
  #
  # root_path <- normalizePath("~/rave_data/raw_dir/testtest/rave-imaging/coregistration/")
  # file.copy(affine_reg_r$fwdtransforms, file.path(root_path, "t1_to_mni_1_0GenericAffine.mat"))
  #
  # fixed_img$to_file(file.path(root_path, "template.nii.gz"))
  # affine_reg_r$warpedmovout$to_file(file.path(root_path, "affine_warpedmovout.nii.gz"))
  # file.copy(syn_reg_r$fwdtransforms[[1]], file.path(root_path, "t1_to_mni_2_deformable_1Warp.nii.gz"))
  # file.copy(syn_reg_r$fwdtransforms[[2]], file.path(root_path, "t1_to_mni_2_deformable_0GenericAffine.mat"))
  #
  # syn_reg_r$warpedmovout$to_file(file.path(root_path, "syn_warpedmovout.nii.gz"))

  re <- structure(
    list(
      description = "Affine+SYN registration to template",
      prefix = outprefix,
      orig = list(
        fixed = "orig_fixed.nii.gz",
        moving = "orig_moving.nii.gz"
      ),
      affine = list(
        warped = "affine_warpedmovout.nii.gz",
        description = "Transform is a 4x4 matrix for both forward and inverse transforms",
        transform = "affine_0GenericAffine.mat"
      ),
      nonlinear = list(
        warped = "deformable_warpedmovout.nii.gz",
        description = "`forward_transforms` include transforms to move from moving to fixed image, and `inverse_transforms` are transforms to move from fixed to moving image.",
        forward_transforms = c(
          "deformable_1Warp.nii.gz",
          "deformable_0GenericAffine.mat"
        ),
        inverse_transforms = c(
          "deformable_0GenericAffine.mat",
          "deformable_1InverseWarp.nii.gz"
        )
      )
    ),
    class = "rpyANTs.halpern_register_template_mri"
  )

  re
}

#' @rdname halpern_preprocess
#' @export
halpern_apply_transform_template_mri <- function(roi_folder, outprefix, verbose = TRUE) {
  # DIPSAUS DEBUG START
  # prefix = "/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/outputs/t1_to_mri_"
  # roi_folder <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/Tower-related/templates 2/Lead-DBS_atlases_MNI_ICBM_2009b_NLIN_ASYM/MNI_ICBM_2009b_NLIN_ASYM/atlases/Functional Connectivity Atlas 7 Networks (Yeo 2011)/mixed"
  # verbose <- TRUE

  prefix <- normalize_path(outprefix, must_work = FALSE)

  if(endsWith(roi_folder, "nii") || endsWith(roi_folder, "nii.gz")) {
    # roi_folder is a file
    rois <- basename(roi_folder)
    roi_folder <- dirname(roi_folder)
  } else {
    rois <- list.files(
      roi_folder,
      pattern = "(nii|nii\\.gz)$",
      full.names = FALSE,
      recursive = TRUE,
      include.dirs = FALSE,
      all.files = FALSE,
      ignore.case = TRUE
    )
  }
  mask_root <- sprintf("%smasks", prefix)
  dir.create(mask_root, showWarnings = FALSE, recursive = TRUE)

  # `warpedmovout`: Moving image warped to space of fixed image.
  # `warpedfixout`: Fixed image warped to space of moving image.
  # `fwdtransforms`: Transforms to move from moving to fixed image.
  # `invtransforms`: Transforms to move from fixed to moving image.
  # gather transforms from template to native
  infiles <- sprintf(
    "%s%s", prefix, c(
      "orig_fixed.nii.gz",
      "orig_moving.nii.gz",
      "affine_0GenericAffine.mat",
      "deformable_0GenericAffine.mat",
      "deformable_1InverseWarp.nii.gz"
    )
  )
  infiles <- normalize_path(infiles, must_work = TRUE)

  native_img <- as_ANTsImage(infiles[[2]])

  # reverse fixed and moving because we want the final image to sit
  # in native space
  transformed <- ants_apply_transforms(
    fixed = native_img,
    moving = infiles[[1]],
    transformlist = infiles[c(3,4,5)],
    interpolator = "nearestNeighbor",
    whichtoinvert = c(TRUE, TRUE, FALSE),
    verbose = verbose
  )

  transformed$to_file(filename = sprintf("%soverall_warpedfixinmov.nii.gz", prefix))

  lapply(rois, function(mask) {
    mask_infile <- normalize_path(file_path(roi_folder, mask), must_work = TRUE)
    mask_outfile <- normalize_path(file_path(mask_root, mask), must_work = FALSE)
    dir.create(dirname(mask_outfile), showWarnings = FALSE, recursive = TRUE)

    mask_out <- ants_apply_transforms(
      fixed = native_img,
      moving = mask_infile,
      transformlist = infiles[c(3,4,5)],
      interpolator = "nearestNeighbor",
      whichtoinvert = c(TRUE, TRUE, FALSE),
      verbose = verbose
    )
    mask_out$to_file(filename = mask_outfile)
    return(NULL)
  })

  invisible()
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
    whichtoinvert <- py_list(as.list(as.logical(whichtoinvert)))
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
