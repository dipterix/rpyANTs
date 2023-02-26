test_that("ANTsTransform generics", {

  # If conda is not set up, skip
  testthat::skip_if_not(rpyANTs:::rpymat_is_setup())
  testthat::skip_if_not(rpyANTs:::ants_available())

  cat("Testing ANTsTransform generics\n")

  cat("Loading numpy\n")
  np <- import("numpy", convert = TRUE)

  cat("Sanity-check on default transform\n")
  y <- ants$new_ants_transform()

  p <- matrix(rnorm(120), nrow = py_to_r(y$dimension))
  vector_trans <- as.matrix(y) %*% rbind(p, 1)
  expected_trans <- rbind(t(np$array(apply(p, 2, y$apply_to_point))), 1)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)
  vector_trans <- as.matrix(y) %*% rbind(p, 0)
  expected_trans <- rbind(t(np$array(apply(p, 2, y$apply_to_vector))), 0)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  cat("Sanity-check on sample data with SyN registration\n")
  ipath1 <- ants$get_ants_data('r16')
  ipath2 <- ants$get_ants_data('r64')

  print(ipath1)
  print(ipath2)
  print(class(ipath1))

  cat("Load images\n")
  fi <- ants$image_read(ants$get_ants_data('r16'))
  mo <- ants$image_read(ants$get_ants_data('r64'))

  print(fi)
  print(dim(fi[]))
  print(mo)
  print(dim(mo[]))

  # Somehow this cannot run on windows
  # # resample to speed up this example
  cat("Resample images for speed\n")
  fi <- ants$resample_image(fi, list(60L,60L), TRUE, 0L)
  mo <- ants$resample_image(mo, list(60L,60L), TRUE, 0L)

  # SDR transform
  cat("Non-linear transform\n")
  transform <- ants_registration(
    fixed=fi, moving=mo, type_of_transform = 'SyN', verbose = TRUE)

  tmp_files <- unique(unlist(c(
    py_to_r(transform$fwdtransforms),
    py_to_r(transform$invtransforms)
  )))
  cat("Transform generated. Needs to remove later\n")
  cat(tmp_files, sep = "\n")
  cat("\n")


  # AffineTransform
  cat("Testing with saved transformation\n")
  y <- ants$read_transform(transform$fwdtransforms[1L])
  z1 <- as_ANTsTransform(y)
  z2 <- as_ANTsTransform(y[], dimension = y$dimension)

  p <- matrix(rnorm(120), nrow = py_to_r(y$dimension))
  vector_trans <- as.matrix(y) %*% rbind(p[,,drop = FALSE], 1)
  expected_trans <- rbind(t(np$array(apply(p, 2, y$apply_to_point))), 1)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  expected_trans <- rbind(t(np$array(apply(p, 2, z1$apply_to_point))), 1)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  expected_trans <- rbind(t(np$array(apply(p, 2, z2$apply_to_point))), 1)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  vector_trans <- as.matrix(y) %*% rbind(p, 0)
  expected_trans <- rbind(t(np$array(apply(p, 2, y$apply_to_vector))), 0)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  expected_trans <- rbind(t(np$array(apply(p, 2, z1$apply_to_vector))), 0)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  expected_trans <- rbind(t(np$array(apply(p, 2, z2$apply_to_vector))), 0)
  expect_equal(vector_trans, expected_trans, tolerance = 1e-4)

  for(f in tmp_files) {
    if(file.exists(f)){ unlink(f) }
  }

})

