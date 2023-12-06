template_urls <- list(
  `mni_icbm152_nlin_asym_09a` = list(
    coord_sys = "MNI152",
    url = "https://github.com/dipterix/threeBrain-sample/releases/download/1.0.1/mni_icbm152_nlin_asym_09a.zip"
  ),
  `mni_icbm152_nlin_asym_09b` = list(
    coord_sys = "MNI152",
    "https://github.com/dipterix/threeBrain-sample/releases/download/1.0.1/mni_icbm152_nlin_asym_09b.zip"
  )
)

BUILTIN_TEMPLATES <- names(template_urls)

ensure_template <- function(name = BUILTIN_TEMPLATES) {
  name <- match.arg(name)
  template_path <- file.path(R_user_dir(package = "rpyANTs", which = "data"), "templates", name)
  if(!dir.exists(template_path)) {
    url <- template_urls[[name]]$url
    f <- tempfile(fileext = ".zip")
    # options("timeout" = 3600)
    utils::download.file(url, destfile = f)
    utils::unzip(f, exdir = dirname(template_path), overwrite = TRUE)
  }
  normalize_path(template_path, must_work = TRUE)
}


#' Process 'T1' image
#' @description
#' Process 'MRI' and align with template brains
#' @param t1_path path to a 'T1' image
#' @param templates template to use; default is \code{'mni_icbm152_nlin_asym_09a'},
#' @param work_path working path, must be a directory
#' @param verbose whether to verbose the progress
#' @returns Nothing will be returned. Please check \code{work_path} for results.
#' @export
t1_preprocess <- function(t1_path, templates = "mni_icbm152_nlin_asym_09a", work_path = ".", verbose = TRUE) {

  templates <- templates[templates %in% BUILTIN_TEMPLATES]
  if(!length(templates)) {
    templates <- BUILTIN_TEMPLATES[[1]]
    message("No valid template is found/chosen, using: ", templates)
  }

  # DIPSAUS DEBUG START
  # t1_path <- "~/rave_data/raw_dir/yael_demo_001/rave-imaging/derivative/MRI_reference.nii.gz"
  # template <- "mni_icbm152_nlin_asym_09a"
  # work_path <- "~/rave_data/raw_dir/Liming/rave-imaging/"
  # verbose <- TRUE
  # templates <- "mni_icbm152_nlin_asym_09a"

  # ---- Preparation -----------------------------------------------------------
  work_path <- dir_create2(work_path)

  # For self-contained RAVE imaging preprocessing, copy the MRI and CT to rave-imaging path
  t1_path2 <- file.path(dir_create2(file.path(work_path, "inputs", "MRI")), "MRI_RAW.nii.gz")
  t1 <- RNifti::readNifti(t1_path, internal = TRUE)
  RNifti::writeNifti(t1, file = t1_path2, compression = 0)

  # ---- Template mapping ------------------------------------------------------
  template_mapping_path <- dir_create2(file.path(work_path, "template_mapping"))
  template_filename <- "T1.nii.gz"
  template_maskname <- "T1_brainmask.nii.gz"

  template_mapping_results <- lapply(templates, function(template_name) {
    template_path <- ensure_template(template_name)
    template_coord_sys <- template_urls[[template_name]]$coord_sys
    template_mapping_prefix <- file_path(template_mapping_path, template_name, sprintf("%s_", template_name))

    # Determine if MNI mask file is absent
    template_mask_path <- file.path(template_path, template_maskname)
    if(!file.exists(template_mask_path)) {
      template_mask_path <- NULL
    }

    # Affine+SYN registration
    message("Registering T1 to template `", template_name, "` (Affine + SYN)")
    re <- halpern_register_template_mri(
      fixed = file.path(template_path, template_filename),
      moving = t1_path2,
      mask = template_mask_path,
      outprefix = template_mapping_prefix,
      verbose = verbose
    )
    re$template <- template_name
    re$template_coord_sys <- template_coord_sys

    # Apply transforms to ROI from template to T1
    if(length(template_mask_path)) {
      halpern_apply_transform_template_mri(
        roi_folder = template_mask_path,
        outprefix = template_mapping_prefix,
        verbose = verbose
      )

      t1_image <- as_ANTsImage(sprintf("%sorig_moving.nii.gz", template_mapping_prefix))
      mask_image <- as_ANTsImage(file_path(sprintf("%smasks", template_mapping_prefix), template_maskname))
      t1_image <- t1_image * mask_image
      t1_image$to_file(sprintf("%sorig_moving_skullstrip.nii.gz", template_mapping_prefix))
      mask_image$to_file(sprintf("%sorig_fixed_mask.nii.gz", template_mapping_prefix))
      re$orig$fixed_mask <- "orig_fixed_mask.nii.gz"
    }

    json <- rpymat::import("json")
    json_string <- to_r(json$dumps(re))
    writeLines(json_string, file.path(template_mapping_path, template_name, "mapping_configuration.json"))

    return(re)
  })


  # ---- Prepare ants folder in case freesurfer is absent (as a backup) ----------
  t1_mri_processing_path <- dir_create2(file.path(work_path, "ants", "mri"))
  t1_transforms_path <- dir_create2(file.path(t1_mri_processing_path, "transforms"))
  first_mapping <- template_mapping_results[[1]]
  template_mapping_prefix <- first_mapping$prefix

  # TODO: change T1 sform so the space is MNI305?
  file.copy(
    from = sprintf("%s%s", first_mapping$prefix, first_mapping$orig$moving),
    to = file.path(t1_mri_processing_path, "T1.nii.gz"),
    overwrite = TRUE
  )

  skullstrip <- sprintf("%sorig_moving_skullstrip.nii.gz", template_mapping_prefix)
  if(file.exists(skullstrip)) {
    file.copy(
      from = skullstrip,
      to = file.path(t1_mri_processing_path, "brain.nii.gz"),
      overwrite = TRUE
    )
  }

  bran_mask <- sprintf("%sorig_fixed_mask.nii.gz", template_mapping_prefix)
  if(file.exists(bran_mask)) {
    file.copy(
      from = bran_mask,
      to = file.path(t1_mri_processing_path, "brainmask.nii.gz"),
      overwrite = TRUE
    )
  }

  # generate t1 to MNI305 transform (RAS)
  template_to_mni305 <- diag(rep(1, 4))
  if( first_mapping$template_coord_sys == "MNI152" ) {
    template_to_mni305 <- solve(MNI305_to_MNI152)
  }
  # get affine transform (TODO: FIXME?)
  t1_to_mni_lps1 <- rpyANTs::as_ANTsTransform(sprintf(
    "%s%s",
    first_mapping$prefix,
    first_mapping$affine$transform
  ))
  t1_to_mni_lps1 <- solve(as.matrix(t1_to_mni_lps1))
  t1_to_mni_lps2 <- rpyANTs::as_ANTsTransform(sprintf(
    "%s%s",
    first_mapping$prefix,
    first_mapping$nonlinear$forward_transforms[[2]]
  ))
  t1_to_mni_lps2 <- solve(as.matrix(t1_to_mni_lps2))

  lps_to_ras <- diag(c(-1, -1, 1, 1))
  t1_to_mni305 <- template_to_mni305 %*% lps_to_ras %*% t1_to_mni_lps2 %*% t1_to_mni_lps1 %*% solve(lps_to_ras)
  t1_to_mni305 <- template_to_mni305 %*% lps_to_ras %*% t1_to_mni_lps1 %*% solve(lps_to_ras)
  # t1_to_mni152 <- MNI305_to_MNI152 %*% t1_to_mni305
  #
  # # calculate rigid transform to MNI152
  # mni152_to_t1 <- solve(t1_to_mni152)
  # for(ii in c(1,2,3)) {
  #   mni152_to_t1[,ii] <- mni152_to_t1[,ii] / norm(mni152_to_t1[,ii], type = "2")
  # }
  # t1_to_mni152_rigid <- solve(mni152_to_t1)

  # t1 <- RNifti::readNifti(t1_path2, internal = TRUE)
  # sform <- RNifti::xform(t1, useQuaternionFirst = TRUE)
  # sform_152 <- t1_to_mni152_rigid %*% sform
  # RNifti::sform(t1) <- sform_152
  # # RNifti::qform(t1) <- sform_152
  # # t1_header <- RNifti::niftiHeader(t1)
  # # t1_header$sform_code <- 1L
  # RNifti::writeNifti(
  #   image = t1,
  #   file = file.path(t1_mri_processing_path, "brain_acpc.nii.gz"),
  #   compression = 0
  # )

  writeLines(
    con = file.path(t1_transforms_path, "talairach.xfm"),
    sprintf(r"(MNI Transform File
%% avi2talxfm

Transform_Type = Linear;
Linear_Transform =
  %s
%s
%s;
)",
  paste(sprintf("%.6f", t1_to_mni305[1, ]), collapse = " "),
  paste(sprintf("%.6f", t1_to_mni305[2, ]), collapse = " "),
  paste(sprintf("%.6f", t1_to_mni305[3, ]), collapse = " ")
    )
  )

  return(invisible())
}


# ct_path <- "~/rave_data/raw_dir/testtest/rave-imaging/coregistration/CT_RAW.nii.gz"
# template_root <- "~/Dropbox (PennNeurosurgery)/RAVE/Samples/Tower-related/templates 2/Lead-DBS_atlases_MNI_ICBM_2009b_NLIN_ASYM/MNI_ICBM_2009b_NLIN_ASYM/"
# template_name <- "MNI_ICBM_2009b_NLIN_ASYM"
# template_coord_sys <- "MNI152"
# template_filename <- "t1_brain.nii.gz"
# template_maskname <- "t1_brain_mask.nii.gz"
# template_atlas <- "atlases/Cortical Functional Networks (Yeo 2011)/"
# freesurfer_home <- raveio::cmd_freesurfer_home()
# freesurfer_skip_reconall <- FALSE
#
# subject_code <- "Liming2"
#
# # ---- Preparation -------------------------------------------------------------
#
# # Create RAVE subject, initialize paths
# subject <- raveio::RAVESubject$new(project_name = "YAEL", subject_code = subject_code, strict = FALSE)
# subject$initialize_paths()
# imaging_root <- subject$imaging_path
#
# # T1 MRI preprocessing
# t1_processing_path <- raveio::dir_create2(file.path(imaging_root, "ants"))
# template_mapping_prefix <- file.path(t1_processing_path, "templates", template_name, sprintf("%s_", template_name))
#
# # ---- Map T1 to MNI template --------------------------------------------------
#
# # For self-contained RAVE imaging preprocessing, copy the MRI and CT to rave-imaging path
# mri_input_dirpath <- raveio::dir_create2(file.path(imaging_root, "inputs", "MRI"))
# file.copy(t1_path, to = file.path(mri_input_dirpath, basename(t1_path)))
#
# ct_input_dirpath <- raveio::dir_create2(file.path(imaging_root, "inputs", "CT"))
# file.copy(ct_path, to = file.path(ct_input_dirpath, basename(ct_path)))
#
# # Determine if MNI mask file is absent
# if(length(template_maskname)) {
#   template_mask_path <- file.path(template_root, template_maskname)
# } else {
#   template_mask_path <- NULL
# }
#
# # Affine+SYN registration
# rpyANTs::halpern_register_template_mri(
#   fixed = file.path(template_root, template_filename),
#   moving = t1_path,
#   mask = template_mask_path,
#   outprefix = template_mapping_prefix,
#   verbose = TRUE
# )
#
# # Apply transforms to ROI from template to T1
# rpyANTs::halpern_apply_transform_template_mri(
#   roi_folder = file.path(template_root, template_atlas),
#   outprefix = template_mapping_prefix,
#   verbose = TRUE
# )
#
# # ---- Prepare ants folder in case freesurfer is absent (as a backup) ----------
# t1_mri_processing_path <- raveio::dir_create2(file.path(t1_processing_path, "mri"))
# t1_transforms_path <- raveio::dir_create2(file.path(t1_processing_path, "mri", "transforms"))
# file.copy(
#   from = sprintf("%sorig_moving.nii.gz", template_mapping_prefix),
#   to = file.path(t1_mri_processing_path, "brain.nii.gz"),
#   overwrite = TRUE
# )
#
# # generate t1 to MNI305 transform (RAS)
# template_to_mni305 <- diag(rep(1, 4))
# if( template_coord_sys == "MNI152" ) {
#   template_to_mni305 <- solve(raveio::MNI305_to_MNI152)
# }
# # get affine transform (TODO: FIXME?)
# t1_to_mni_lps <- rpyANTs::as_ANTsTransform(sprintf("%saffine_0GenericAffine.mat", template_mapping_prefix))
# t1_to_mni_lps <- solve(as.matrix(t1_to_mni_lps))
# lps_to_ras <- diag(c(-1, -1, 1, 1))
# t1_to_mni305 <- template_to_mni305 %*% lps_to_ras %*% t1_to_mni_lps %*% solve(lps_to_ras)
#
# writeLines(
#   con = file.path(t1_transforms_path, "talairach.xfm"),
#   sprintf(r"(MNI Transform File
# %% avi2talxfm
#
# Transform_Type = Linear;
# Linear_Transform =
#   %s
# %s
# %s;
# )",
# paste(sprintf("%.6f", t1_to_mni305[1, ]), collapse = " "),
# paste(sprintf("%.6f", t1_to_mni305[2, ]), collapse = " "),
# paste(sprintf("%.6f", t1_to_mni305[3, ]), collapse = " ")
#   )
# )
#
# # Apply the mask to t1, generating "brain.finalsurfs.nii.gz"
# if(length(template_mask_path)) {
#   rpyANTs::halpern_apply_transform_template_mri(
#     roi_folder = template_mask_path,
#     outprefix = template_mapping_prefix,
#     verbose = TRUE
#   )
#   native_mask <- sprintf("%smasks/%s", template_mapping_prefix, basename(template_mask_path))
#   file.copy(native_mask, file.path(t1_mri_processing_path, "brainmask.nii.gz"))
#
#   t1_img <- rpyANTs::as_ANTsImage(file.path(t1_mri_processing_path, "brain.nii.gz"))
#   mask_img <- rpyANTs::as_ANTsImage(file.path(t1_mri_processing_path, "brainmask.nii.gz"))
#   t1_img <- t1_img * mask_img
#   t1_img$to_file(rpyANTs:::file_path(t1_mri_processing_path, "brain.finalsurfs.nii.gz"))
# }
#
#
# # ---- Apply CT-MRI coregistration ---------------------------------------------
# if(length(ct_path) == 1) {
#   coreg_path <- raveio::dir_create2(file.path(imaging_root, "coregistration"))
#   rpyANTs::halpern_register_ct_mri(
#     fixed = ct_path,
#     moving = t1_path,
#     outprefix = file.path(coreg_path, "t1_to_ct_"),
#     fixed_is_ct = TRUE,
#     verbose = TRUE
#   )
# }
#
# # ---- Prepare freesurfer folder (autorecon1) ----------------------------------
# fs_path <- file.path(imaging_root, "fs")
# fs_path_isvalid <- threeBrain::check_freesurfer_path(fs_path, check_volume = TRUE, check_surface = FALSE)
# if(isFALSE(fs_path_isvalid)) {
#   # Need to run FreeSurfer
#   raveio::cmd_run_recon_all(
#     subject = subject,
#     mri_path = t1_path,
#     args = "-autorecon1",
#     work_path = NULL,
#     overwrite = TRUE,
#     command_path = freesurfer_home,
#     dry_run = FALSE,
#     verbose = TRUE
#   )
#
#   slice_file <- file.path(t1_mri_processing_path,
#                           c("brain.finalsurfs.nii.gz", "brain.nii.gz"))
#   slice_file <- slice_file[file.exists(slice_file)][[1]]
#   file.copy(slice_file, file.path(fs_path, "mri", "rave_slices.nii.gz"))
#
# }
#
# # ---- Finish up the rest of recon-all in the background -----------------------
# if(!freesurfer_skip_reconall) {
#   # Finish the rest of recon-all
#   dipsaus::rs_exec({
#     raveio::cmd_run_recon_all(
#       subject = subject,
#       mri_path = t1_path,
#       args = "-all",
#       work_path = NULL,
#       overwrite = FALSE,
#       command_path = freesurfer_home,
#       dry_run = FALSE,
#       verbose = TRUE
#     )
#   }, name = sprintf("recon-all (%s)", subject$subject_code),
#   focus_on_console = TRUE)
# }
#
# # ---- Launch RAVE -------------------------------------------------------------
# message("Running FreeSurfer recon-all in the background. Launching localization... Please do NOT close RStudio. Wait for recon-all to finish.")
# rave::start_yael(as_job = TRUE)
#
#
# # ---- Apply transforms to electrodes
# # PLEASE RUN THIS AFTER ELECTRODE LOCALIZATION!!!
#
# brain <- raveio::rave_brain(subject)
# electrode_table <- brain$electrodes$raw_table
# t1_ras <- electrode_table[, c("T1R", "T1A", "T1S")]
# invalids <- rowSums(t1_ras^2) == 0
#
# transforms <- rpyANTs:::normalize_path(
#   sprintf("%s%s", template_mapping_prefix, c("affine_0GenericAffine.mat",
#                                              "deformable_0GenericAffine.mat",
#                                              "deformable_1InverseWarp.nii.gz"))
# )
#
# # ANTs uses LPS instead of RAS
# t1_lps <- data.frame(
#   x = -t1_ras$T1R,
#   y = -t1_ras$T1A,
#   z = t1_ras$T1S
# )
#
# mni_lps <- rpyANTs::ants_apply_transforms_to_points(
#   dim = 3,
#   points = t1_lps,
#   transformlist = transforms,
#   whichtoinvert = c(TRUE, TRUE, FALSE),
#   verbose = FALSE
# )
# mni_lps <- data.matrix(mni_lps)
# mni_lps[invalids, ] <- 0
#
# if(template_coord_sys == "MNI152") {
#   electrode_table$MNI152_x <- -mni_lps[,1]
#   electrode_table$MNI152_y <- -mni_lps[,2]
#   electrode_table$MNI152_z <- mni_lps[,3]
#
#   mni305 <- solve(raveio::MNI305_to_MNI152) %*% rbind(-mni_lps[,1], -mni_lps[,2], mni_lps[,3], 1)
#   electrode_table$MNI305_x <- mni305[1, ]
#   electrode_table$MNI305_y <- mni305[2, ]
#   electrode_table$MNI305_z <- mni305[3, ]
# } else {
#   electrode_table$MNI305_x <- -mni_lps[,1]
#   electrode_table$MNI305_y <- -mni_lps[,2]
#   electrode_table$MNI305_z <- mni_lps[,3]
#
#   mni152 <- raveio::MNI305_to_MNI152 %*% rbind(-mni_lps[,1], -mni_lps[,2], mni_lps[,3], 1)
#   electrode_table$MNI152_x <- mni152[1, ]
#   electrode_table$MNI152_y <- mni152[2, ]
#   electrode_table$MNI152_z <- mni152[3, ]
# }
# brain$set_electrodes(electrode_table)
#
# raveio::save_meta2(data = electrode_table, meta_type = "electrodes", project_name = subject$project_name, subject_code = subject$subject_code)
#
#
# # threeBrain::merge_brain(brain)$plot()
#
# # get atlas
# brain <- raveio::rave_brain(subject)
#
#
#
#
#
#
#
