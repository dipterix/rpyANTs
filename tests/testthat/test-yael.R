library(testthat)

test_that("YAEL session labels follow the image types", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  rpyants <- rpyANTs:::load_rpyants()

  # Only file names are formatted here, nothing is written to `work_path`
  work_path <- tempfile()
  on.exit({ unlink(work_path, recursive = TRUE) }, add = TRUE)

  preop_types <- c("T1w", "T2w", "FLAIR", "preopCT")
  postop_types <- c("CT", "postopCT", "postopT1w", "postopT2w", "postopFLAIR")

  yael <- rpyants$registration$YAELPreprocess(
    "Demo", work_path, as.list(c(preop_types, postop_types)))

  input_name <- function(type, ...) {
    path <- yael$format_path(folder = "inputs/anat", name = type,
                             ext = "nii.gz", relative = TRUE, ...)
    basename(py_to_r(path))
  }

  for (type in preop_types) {
    expect_equal(
      input_name(type),
      sprintf("sub-Demo_ses-preop_desc-preproc_%s.nii.gz", type)
    )
  }

  for (type in postop_types) {
    expect_equal(
      input_name(type),
      sprintf("sub-Demo_ses-postop_desc-preproc_%s.nii.gz", type)
    )
  }

  # The session is only inferred when it is not given
  expect_equal(
    input_name("postopT1w", ses = "implant"),
    "sub-Demo_ses-implant_desc-preproc_postopT1w.nii.gz"
  )

})

test_that("YAEL finds `postop` images saved with the previous session label", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  rpyants <- rpyANTs:::load_rpyants()

  work_path <- tempfile()
  on.exit({ unlink(work_path, recursive = TRUE) }, add = TRUE)

  # rpyANTs <= 0.0.6 labeled `postop*` images as `ses-preop`; the files are
  # empty since only their names are parsed
  input_dir <- file.path(work_path, "inputs", "anat")
  aligned_dir <- file.path(work_path, "coregistration", "anat")
  dir.create(input_dir, recursive = TRUE)
  dir.create(aligned_dir, recursive = TRUE)

  old_input <- "sub-Demo_ses-preop_desc-preproc_postopT1w.nii.gz"
  old_aligned <- "sub-Demo_ses-preop_space-scanner_desc-preproc_postopT1w.nii.gz"
  new_aligned <- "sub-Demo_ses-postop_space-scanner_desc-preproc_postopT1w.nii.gz"
  file.create(file.path(input_dir, old_input))
  file.create(file.path(aligned_dir, old_aligned))

  yael <- rpyants$registration$YAELPreprocess(
    "Demo", work_path, list("T1w", "postopT1w"))

  aligned_name <- function() {
    mapping <- py_to_r(yael$get_native_mapping("postopT1w", relative = TRUE))
    basename(mapping$mappings$postopT1w_in_T1w)
  }

  expect_equal(basename(py_to_r(yael$input_image_path("postopT1w"))), old_input)
  expect_equal(aligned_name(), old_aligned)

  # Registering the image again writes the aligned image with the new label,
  # next to the old one: the mapping must use the new one, also when the old
  # one is listed last (as in sorted order)
  file.create(file.path(aligned_dir, new_aligned))
  expect_equal(with_sorted_listdir(aligned_name()), new_aligned)

})

test_that("YAEL returns the aligned image as it is named on disk", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  rpyants <- rpyANTs:::load_rpyants()

  work_path <- tempfile()
  on.exit({ unlink(work_path, recursive = TRUE) }, add = TRUE)

  # The files spell the subject code in lower case. On file systems that ignore
  # the case of file names, a file named with the given subject code "exists"
  input_dir <- file.path(work_path, "inputs", "anat")
  aligned_dir <- file.path(work_path, "coregistration", "anat")
  dir.create(input_dir, recursive = TRUE)
  dir.create(aligned_dir, recursive = TRUE)

  aligned <- "sub-demo_ses-postop_space-scanner_desc-preproc_postopT1w.nii.gz"
  file.create(file.path(input_dir, "sub-demo_ses-postop_desc-preproc_postopT1w.nii.gz"))
  file.create(file.path(aligned_dir, aligned))

  yael <- rpyants$registration$YAELPreprocess(
    "Demo", work_path, list("T1w", "postopT1w"))

  mapping <- py_to_r(yael$get_native_mapping("postopT1w", relative = TRUE))
  expect_equal(basename(mapping$mappings$postopT1w_in_T1w), aligned)

})
