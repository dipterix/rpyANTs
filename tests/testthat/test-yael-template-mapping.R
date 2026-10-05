library(testthat)

test_that("YAEL template mapping handles stray files and renamed subjects", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  rpyants <- rpyANTs:::load_rpyants()

  work_path <- tempfile()
  on.exit({ unlink(work_path, recursive = TRUE) }, add = TRUE)
  transform_dir <- file.path(work_path, "normalization", "transformations")
  input_dir <- file.path(work_path, "inputs", "anat")
  dir.create(transform_dir, recursive = TRUE)
  dir.create(input_dir, recursive = TRUE)

  # only the file names are parsed, the files can be empty
  transform_files <- function(sub, template = "MNI152NLin2009bAsym") {
    c(
      sprintf("sub-%s_from-T1w_to-%s_desc-affine+SyN_ants0.mat", sub, template),
      sprintf("sub-%s_from-T1w_to-%s_desc-affine+SyN_ants1.nii.gz", sub, template),
      sprintf("sub-%s_from-%s_to-T1w_desc-affine+SyN_ants0.nii.gz", sub, template),
      sprintf("sub-%s_from-%s_to-T1w_desc-affine+SyN_ants1.mat", sub, template)
    )
  }
  set_native_image <- function(sub) {
    unlink(list.files(input_dir, full.names = TRUE))
    file.create(file.path(input_dir, sprintf("sub-%s_ses-preop_desc-preproc_T1w.nii.gz", sub)))
  }
  get_mapping <- function(subject_code) {
    yael <- rpyants$registration$YAELPreprocess(subject_code, work_path)
    py_to_r(yael$get_template_mapping(
      template_name = "MNI152NLin2009bAsym", relative = TRUE))
  }
  transform_names <- function(mapping, which) {
    basename(unlist(mapping[[which]]$transformlist))
  }

  # files without `from`/`to` entities are ignored
  file.create(file.path(transform_dir, c(transform_files("Demo"),
                                         "sub-Demo_desc-foo_ants0.mat")))
  set_native_image("Demo")
  mapping <- get_mapping("Demo")
  expect_equal(transform_names(mapping, "native_to_template"),
               transform_files("Demo")[c(2, 1)])
  expect_equal(transform_names(mapping, "template_to_native"),
               transform_files("Demo")[c(4, 3)])

  # the subject code is matched case-insensitively (without the native image,
  # so the renamed-folder rule below cannot be the reason)
  set_native_image(NULL)
  expect_equal(transform_names(get_mapping("demo"), "native_to_template"),
               transform_files("Demo")[c(2, 1)])
  set_native_image("Demo")

  # renamed subject folder: the transforms and the native image share
  # the old subject code
  expect_equal(transform_names(get_mapping("Renamed"), "template_to_native"),
               transform_files("Demo")[c(4, 3)])

  # ... but not when the native image was imported under the new code, or
  # when it is missing
  set_native_image("Renamed")
  expect_null(get_mapping("Renamed"))
  set_native_image(NULL)
  expect_null(get_mapping("Renamed"))

  # transforms from several other subjects are ambiguous
  set_native_image("Demo")
  file.create(file.path(transform_dir, transform_files("Other")))
  expect_null(get_mapping("Renamed"))

  # the subject's own transforms are preferred
  expect_equal(transform_names(get_mapping("Other"), "native_to_template"),
               transform_files("Other")[c(2, 1)])

  # on case-sensitive file systems, two spellings of the subject code are
  # ambiguous rather than merged
  probe <- file.path(work_path, "case-probe")
  file.create(probe)
  if (!file.exists(toupper(probe))) {
    file.create(file.path(transform_dir, transform_files("DEMO")))
    expect_null(get_mapping("demo"))
    expect_equal(transform_names(get_mapping("Demo"), "native_to_template"),
                 transform_files("Demo")[c(2, 1)])
  }
})
