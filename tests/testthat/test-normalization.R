library(testthat)

test_that("Stage cache keeps both warped images of a registration", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()
  StageContext <- rpyants$utils$cache$StageContext

  root <- tempfile()
  dir.create(root)
  on.exit({ unlink(root, recursive = TRUE) }, add = TRUE)

  # Toy registration result: constant images, and placeholder transform files
  # (the cache moves the transform files without reading them)
  transforms <- file.path(root, c(
    "toy1Warp.nii.gz", "toy0GenericAffine.mat", "toy1InverseWarp.nii.gz"))
  toy_result <- function(...) {
    for (f in transforms) { writeLines("toy", f) }
    images <- list(...)
    py_dict(
      keys = as.list(c(names(images), "fwdtransforms", "invtransforms")),
      values = c(
        unname(images),
        list(as.list(transforms[c(1, 2)]), as.list(transforms[c(2, 3)]))
      ),
      convert = FALSE
    )
  }

  # Python: `with StageContext(...) as (ctx, cached): ctx.result = result`
  # Returns `cached`, and stores `result` when it is given
  enter_stage <- function(result = NULL) {
    ctx <- StageContext("toy", "registration", root, verbose = FALSE)
    cached <- py_get_item(ctx$`__enter__`(), 1L)
    if (!is.null(result)) {
      py_set_attr(ctx, "result", result)
    }
    ctx$`__exit__`(NULL, NULL, NULL)
    cached
  }

  cached <- enter_stage(toy_result(
    warpedmovout = ants$make_image(c(4L, 4L, 4L), voxval = 1),
    warpedfixout = ants$make_image(c(4L, 4L, 4L), voxval = 2)
  ))
  expect_null(py_to_r(cached))

  # The second time the stage is entered, the result comes from the cache
  cached <- enter_stage()

  expect_setequal(
    names(py_to_r(cached)),
    c("warpedmovout", "warpedfixout", "fwdtransforms", "invtransforms")
  )
  expect_equal(py_to_r(py_get_item(cached, "warpedmovout")$sum()), 64)
  expect_equal(py_to_r(py_get_item(cached, "warpedfixout")$sum()), 128)
  expect_equal(
    basename(unlist(py_to_r(py_get_item(cached, "fwdtransforms")))),
    c("toy_fwdtransforms_ants0.nii.gz", "toy_fwdtransforms_ants1.mat")
  )
  expect_equal(
    basename(unlist(py_to_r(py_get_item(cached, "invtransforms")))),
    c("toy_invtransforms_ants0.mat", "toy_invtransforms_ants1.nii.gz")
  )

  # A registration stored again without `warpedfixout` does not get the image
  # of the previous one (removing a cached file makes the stage run again)
  unlink(file.path(root, "toy_warped_mov.nii.gz"))
  cached <- enter_stage(toy_result(
    warpedmovout = ants$make_image(c(4L, 4L, 4L), voxval = 3)
  ))
  expect_null(py_to_r(cached))

  cached <- enter_stage()
  expect_setequal(
    names(py_to_r(cached)),
    c("warpedmovout", "fwdtransforms", "invtransforms")
  )
  expect_equal(py_to_r(py_get_item(cached, "warpedmovout")$sum()), 192)

})

# Settings of the final registration in `normalization_with_atropos`
syn_aggro_settings <- function() {
  list(
    grad_step = 0.15, flow_sigma = 3.5, total_sigma = 0,
    aff_metric = "mattes", aff_sampling = 32L, aff_random_sampling_rate = 0.2,
    syn_metric = "CC", syn_sampling = 4L,
    reg_iterations = tuple(100L, 70L, 50L, 0L),
    verbose = FALSE
  )
}

# Reads the displacement field (in pixels) of a registration and removes the
# transform files from the temporary directory
read_warp <- function(res) {
  fwd <- unlist(py_to_r(py_get_item(res, "fwdtransforms")))
  inv <- unlist(py_to_r(py_get_item(res, "invtransforms")))
  on.exit({ unlink(c(fwd, inv)) })
  py_to_r(ants$image_read(fwd[[1]])$numpy())
}

test_that("`registration_syn_aggro` uses the additional metrics", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()
  registration_syn_aggro <- rpyants$registration$normalization$registration_syn_aggro

  # The main images are identical, so they give nothing to deform. Only the
  # additional channel differs: a blob sits at different places
  fixed <- toy_disks(c(48, 48, 30, 100), c(48, 48, 16, 60))
  moving <- fixed$clone()
  extra_fixed <- toy_disks(c(36, 48, 7, 1))
  extra_moving <- toy_disks(c(60, 48, 7, 1))

  res <- do.call(registration_syn_aggro, c(
    list(fixed = fixed, moving = moving),
    syn_aggro_settings()
  ))
  expect_lt(max(abs(read_warp(res))), 1)

  res <- do.call(registration_syn_aggro, c(
    list(
      fixed = fixed, moving = moving,
      multivariate_extras = list(
        list("MI", extra_fixed, extra_moving, 0.5, "32,Random,0.25")
      )
    ),
    syn_aggro_settings()
  ))

  # Same layout as 'SyNAggro': [warp, affine] and [affine, inverse warp]
  fwd <- unlist(py_to_r(py_get_item(res, "fwdtransforms")))
  inv <- unlist(py_to_r(py_get_item(res, "invtransforms")))
  expect_length(fwd, 2)
  expect_match(fwd[[1]], "1Warp\\.nii\\.gz$")
  expect_match(fwd[[2]], "0GenericAffine\\.mat$")
  expect_length(inv, 2)
  expect_match(inv[[1]], "0GenericAffine\\.mat$")
  expect_match(inv[[2]], "1InverseWarp\\.nii\\.gz$")
  expect_true(all(file.exists(c(fwd, inv))))

  # The additional channel deforms the image
  expect_gt(max(abs(read_warp(res))), 3)

})

test_that("`registration_syn_aggro` without additional metrics registers like 'SyNAggro'", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()
  registration_syn_aggro <- rpyants$registration$normalization$registration_syn_aggro

  fixed <- toy_disks(c(48, 48, 30, 100), c(48, 48, 16, 60))
  moving <- toy_disks(c(52, 44, 27, 100), c(54, 46, 13, 60))

  # An empty list still runs the affine and SyN stages separately
  res <- do.call(registration_syn_aggro, c(
    list(fixed = fixed, moving = moving, multivariate_extras = list()),
    syn_aggro_settings()
  ))
  warp <- read_warp(res)

  res <- do.call(ants$registration, c(
    list(fixed = fixed, moving = moving, type_of_transform = "SyNAggro"),
    syn_aggro_settings()
  ))
  warp_syn_aggro <- read_warp(res)

  # Displacements are up to 4 pixels; two runs of 'SyNAggro' itself differ by
  # 0.02 to 0.04 pixels on average because its affine stage samples randomly
  expect_gt(max(abs(warp_syn_aggro)), 2)
  expect_lt(mean(abs(warp - warp_syn_aggro)), 0.25)

})

test_that("`registration_syn_aggro` masks its affine stage only with `mask_all_stages`", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()
  registration_syn_aggro <- rpyants$registration$normalization$registration_syn_aggro

  fixed <- toy_disks(c(48, 48, 30, 100), c(48, 48, 16, 60))
  moving <- toy_disks(c(52, 44, 27, 100), c(54, 46, 13, 60))
  mask <- ants$get_mask(fixed)

  spy <- spy_ants_registration()
  on.exit({ spy$restore() }, add = TRUE)

  # Returns the calls of the two stages
  register <- function(mask_all_stages) {
    res <- do.call(registration_syn_aggro, c(
      list(fixed = fixed, moving = moving, multivariate_extras = list(),
           mask = mask, mask_all_stages = mask_all_stages),
      syn_aggro_settings()
    ))
    read_warp(res)
    utils::tail(spy$calls(), 2)
  }

  # 'SyNAggro' always masks its SyN stage
  calls <- register(mask_all_stages = FALSE)
  expect_equal(calls$type, c("Affine", "SyNOnly"))
  expect_equal(calls$masked, c(FALSE, TRUE))

  calls <- register(mask_all_stages = TRUE)
  expect_equal(calls$type, c("Affine", "SyNOnly"))
  expect_equal(calls$masked, c(TRUE, TRUE))

})

test_that("`registration_syn_aggro` keeps the affine settings of 'SyNAggro' and only leaves the transforms of the registration", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()
  registration_syn_aggro <- rpyants$registration$normalization$registration_syn_aggro

  fixed <- toy_disks(c(48, 48, 30, 100), c(48, 48, 16, 60))
  moving <- toy_disks(c(52, 44, 27, 100), c(54, 46, 13, 60))

  # `out_dir` is for the requested outputs; `tmp_dir` shows what 'ANTsPy'
  # writes besides them
  root <- tempfile()
  out_dir <- file.path(root, "out")
  tmp_dir <- file.path(root, "tmp")
  dir.create(out_dir, recursive = TRUE)
  on.exit({ unlink(root, recursive = TRUE) }, add = TRUE)
  restore_tempdir <- use_python_tempdir(tmp_dir)
  on.exit({ restore_tempdir() }, add = TRUE)

  spy <- spy_ants_registration()
  on.exit({ spy$restore() }, add = TRUE)

  # 'SyNAggro' ignores the settings of the other affine transforms
  # (`aff_iterations`, ...): its affine stage keeps its own when they are given
  res <- do.call(registration_syn_aggro, c(
    list(fixed = fixed, moving = moving, multivariate_extras = list(),
         outprefix = file.path(out_dir, "toy_"),
         aff_iterations = tuple(10L, 0L, 0L, 0L),
         restrict_transformation = tuple(1, 1)),
    syn_aggro_settings()
  ))
  calls <- spy$calls()
  expect_equal(calls$type, c("Affine", "SyNOnly"))
  expect_equal(calls$affine[[1]], "2100x1200x1200x100 4x2x2x1 3x2x1x0")

  # `restrict_transformation` only applies to the first stage of a
  # registration, which is the affine stage of 'SyNAggro'
  expect_equal(calls$restricted, c(TRUE, FALSE))

  # The SyN stage starts from the affine, given as a file path: 'ANTsPy' < 0.4
  # does not take a list
  expect_equal(calls$initial, c("NoneType", "str"))

  # Only the transforms of the registration are written, none is missing
  expect_equal(
    basename(unlist(py_to_r(py_get_item(res, "fwdtransforms")))),
    c("toy_1Warp.nii.gz", "toy_0GenericAffine.mat")
  )
  expect_equal(
    basename(unlist(py_to_r(py_get_item(res, "invtransforms")))),
    c("toy_0GenericAffine.mat", "toy_1InverseWarp.nii.gz")
  )
  expect_setequal(
    list.files(out_dir),
    c("toy_0GenericAffine.mat", "toy_1Warp.nii.gz", "toy_1InverseWarp.nii.gz")
  )
  expect_length(list.files(tmp_dir), 0)

  # With composite transforms, `ants.registration` returns one file path per
  # direction instead of lists
  res <- do.call(registration_syn_aggro, c(
    list(fixed = fixed, moving = moving, multivariate_extras = list(),
         outprefix = file.path(out_dir, "composite_"),
         write_composite_transform = TRUE),
    syn_aggro_settings()
  ))
  expect_equal(
    basename(py_to_r(py_get_item(res, "fwdtransforms"))),
    "composite_Composite.h5"
  )
  expect_equal(
    basename(py_to_r(py_get_item(res, "invtransforms"))),
    "composite_InverseComposite.h5"
  )
  expect_setequal(
    list.files(out_dir, pattern = "^composite_"),
    c("composite_Composite.h5", "composite_InverseComposite.h5")
  )
  expect_length(list.files(tmp_dir), 0)

  # When the SyN stage fails (an additional metric needs five items), the
  # affine of the first stage is removed as well
  expect_error(
    do.call(registration_syn_aggro, c(
      list(fixed = fixed, moving = moving,
           multivariate_extras = list(list("MI", fixed))),
      syn_aggro_settings()
    ))
  )
  expect_length(list.files(tmp_dir), 0)

})

test_that("`registration_syn_aggro` takes additional images on another grid", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()
  registration_syn_aggro <- rpyants$registration$normalization$registration_syn_aggro

  # The probability maps of a template can be finer than the template image:
  # here the fixed blob has twice the resolution of the other images
  fixed <- toy_disks(c(48, 48, 30, 100), c(48, 48, 16, 60))
  moving <- fixed$clone()
  extra_fixed <- toy_disks(c(36, 48, 7, 1), spacing = 0.5)
  extra_moving <- toy_disks(c(60, 48, 7, 1))
  expect_equal(unlist(py_to_r(extra_fixed$shape)), c(192, 192))

  res <- do.call(registration_syn_aggro, c(
    list(
      fixed = fixed, moving = moving,
      multivariate_extras = list(
        list("MI", extra_fixed, extra_moving, 0.5, "32,Random,0.25")
      )
    ),
    syn_aggro_settings()
  ))
  warp <- read_warp(res)

  # The displacement field is on the grid of `fixed`, and the additional
  # channel still deforms the image (less than 0.5 pixels without it)
  expect_equal(dim(warp), c(96, 96, 2))
  expect_gt(max(abs(warp)), 1.5)

})

test_that("Normalization needs as many template images as subject images, or one", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  rpyants <- rpyANTs:::load_rpyants()

  # The images are paired to make the additional metrics. The paths do not
  # exist: the lengths are checked before any image is read
  expect_error(
    rpyants$registration$normalization$normalization_with_atropos(
      fix_path = list("template/T1.nii.gz", "template/T2.nii.gz"),
      mov_paths = list("subject_T1.nii.gz"),
      working_path = tempfile(), use_antspynet = FALSE, verbose = FALSE
    ),
    "must have the same length"
  )

})

test_that("Normalization without `antspynet` aligns a toy head and caches the final registration", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  # FIXME: No check on Github MacOS due to ImageIO issue
  testthat::skip_if(nzchar(Sys.getenv("GITHUB_OUTPUT")))

  rpyants <- rpyANTs:::load_rpyants()

  root <- tempfile()
  dir.create(root)
  on.exit({ unlink(root, recursive = TRUE) }, add = TRUE)

  toy <- toy_head_dataset(root)
  work_path <- file.path(root, "work")

  # `ants.atropos` does not remove its temporary files: keep them in `root`
  restore_tempdir <- use_python_tempdir(file.path(root, "tmp"))
  on.exit({ restore_tempdir() }, add = TRUE)

  # Two images of the subject, with their weights given as a tuple
  normalize <- function() {
    res <- NULL
    # the stages print their cache status
    reticulate::py_capture_output({
      res <- rpyants$registration$normalization$normalization_with_atropos(
        fix_path = toy$template_t1,
        mov_paths = list(toy$subject_t1, toy$subject_t2),
        working_path = work_path, weights = tuple(1, 0.9),
        use_antspynet = FALSE, verbose = FALSE
      )
    })
    res
  }
  transform_paths <- function(res) {
    list(
      fwd = unlist(py_to_r(py_get_item(res, "fwdtransforms"))),
      inv = unlist(py_to_r(py_get_item(res, "invtransforms")))
    )
  }

  spy <- spy_ants_registration()
  on.exit({ spy$restore() }, add = TRUE)

  res <- normalize()
  paths <- transform_paths(res)

  # After the initial registration, the final one runs the two stages of
  # 'SyNAggro', both with masks: the affine stage with the settings of
  # 'SyNAggro', then the SyN stage with three additional metrics, which are
  # the second image (weight 0.9), and the probability maps of CSF (0.1) and
  # of the deep structures (0.5)
  calls <- spy$calls()
  expect_equal(calls$type, c("SyNabp", "Affine", "SyNOnly"))
  expect_equal(calls$masked, c(FALSE, TRUE, TRUE))
  expect_equal(calls$affine, c("", "2100x1200x1200x100 4x2x2x1 3x2x1x0", ""))
  expect_equal(calls$weights, c("", "", "0.9, 0.1, 0.5"))

  # [warp, affine] and [affine, inverse warp], stored in the working directory
  expect_equal(
    basename(paths$fwd),
    c("SyN_w_atropos_fwdtransforms_ants0.nii.gz",
      "SyN_w_atropos_fwdtransforms_ants1.mat")
  )
  expect_equal(
    basename(paths$inv),
    c("SyN_w_atropos_invtransforms_ants0.mat",
      "SyN_w_atropos_invtransforms_ants1.nii.gz")
  )
  expect_true(all(
    file.exists(file.path(work_path, basename(c(paths$fwd, paths$inv))))
  ))

  # The toy subject is aligned to the toy template: every tissue overlaps more
  template_seg <- ants$image_read(toy$template_seg)
  subject_seg <- ants$image_read(toy$subject_seg)
  seg_before <- ants$resample_image_to_target(
    subject_seg, template_seg, interp_type = "nearestNeighbor")
  seg_after <- ants$apply_transforms(
    fixed = template_seg, moving = subject_seg,
    transformlist = as.list(paths$fwd), interpolator = "nearestNeighbor")
  dice_before <- toy_dice(py_to_r(template_seg$numpy()), py_to_r(seg_before$numpy()))
  dice_after <- toy_dice(py_to_r(template_seg$numpy()), py_to_r(seg_after$numpy()))
  expect_true(all(dice_after > dice_before))
  expect_gt(min(dice_after), 0.6)

  # Normalizing again restores both registrations from the cache instead of
  # running them again, with both warped images
  res <- normalize()
  expect_equal(nrow(spy$calls()), 3)
  expect_equal(transform_paths(res), paths)
  expect_true(all(c("warpedmovout", "warpedfixout") %in% names(py_to_r(res))))

})
