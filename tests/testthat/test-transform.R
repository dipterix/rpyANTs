library(testthat)

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



test_that("Check `get_xform` is consistent with RNifti", {

  mat_to_quaternion <- function (m) {
    m <- m[1:3, 1:3]
    m <- apply(m, 2, function(x) {
      l2 <- sum(x^2)
      if (l2 > 0) {
        x <- x/sqrt(l2)
      }
      x
    })
    m11 <- m[1, 1]
    m21 <- m[2, 1]
    m31 <- m[3, 1]
    m12 <- m[1, 2]
    m22 <- m[2, 2]
    m32 <- m[3, 2]
    m13 <- m[1, 3]
    m23 <- m[2, 3]
    m33 <- m[3, 3]
    trace <- m11 + m22 + m33
    if (trace > 0) {
      s <- 0.5/sqrt(trace + 1)
      w <- 0.25/s
      x <- (m32 - m23) * s
      y <- (m13 - m31) * s
      z <- (m21 - m12) * s
    }
    else if (m11 > m22 && m11 > m33) {
      s <- 2 * sqrt(1 + m11 - m22 - m33)
      w <- (m32 - m23)/s
      x <- 0.25 * s
      y <- (m12 + m21)/s
      z <- (m13 + m31)/s
    }
    else if (m22 > m33) {
      s <- 2 * sqrt(1 + m22 - m11 - m33)
      w <- (m13 - m31)/s
      x <- (m12 + m21)/s
      y <- 0.25 * s
      z <- (m23 + m32)/s
    }
    else {
      s <- 2 * sqrt(1 + m33 - m11 - m22)
      w <- (m21 - m12)/s
      x <- (m13 + m31)/s
      y <- (m23 + m32)/s
      z <- 0.25 * s
    }
    c(x = x, y = y, z = z, w = w)
  }

  io_write_nii.array <- function (x, con, vox2ras = NULL, ...) {
    if (!is.matrix(vox2ras)) {
      warning("`io_write_nii.array`: `vox2ras` is missing, using identity matrix. Please specify voxel-to-RAS transform (4x4 matrix).")
      vox2ras <- diag(1, 4)
    }
    stopifnot(is.matrix(vox2ras) && nrow(vox2ras) == 4 && ncol(vox2ras) == 4)
    quaternion <- mat_to_quaternion(vox2ras)
    shape <- dim(x)
    stopifnot(length(shape) %in% c(3, 4))
    if (length(shape) > 3) {
      nframes <- shape[[4]]
    } else {
      nframes <- 1
    }
    m33 <- vox2ras[1:3, 1:3]
    pixdim <- sqrt(colSums(m33^2))
    pixdim <- c(sign(det(m33)), pixdim, nframes, 0, 0, 0)
    pixdim <- as.double(pixdim)
    if (nframes == 1 && length(shape) != 3) {
      pixdim[[5]] <- 0
      shape <- shape[1:3]
      x <- array(x[seq_len(prod(shape))], dim = shape)
    }
    x[is.na(x)] <- 0
    rg <- range(x)
    if (all(x - round(x) == 0)) {
      if (rg[[1]] >= 0 && rg[[2]] <= 255) {
        datatype_code <- 2L
        bitpix <- 8L
        storage.mode(x) <- "integer"
      } else if (rg[[1]] >= -32768 && rg[[2]] <= 32768) {
        datatype_code <- 4L
        bitpix <- 16L
        storage.mode(x) <- "integer"
      } else if (rg[[1]] >= -2147483648 && rg[[2]] <= 2147483648) {
        datatype_code <- 8L
        bitpix <- 32L
        storage.mode(x) <- "integer"
      } else {
        bitpix <- 32L
        datatype_code <- 16L
      }
    } else {
      bitpix <- 32L
      datatype_code <- 16L
    }
    nii <- RNifti::niftiHeader(x)
    nii$datatype <- datatype_code
    nii$bitpix <- bitpix
    nii$pixdim <- pixdim
    nii$xyzt_units <- 10L
    nii$qform_code <- 1L
    nii$sform_code <- 1L
    qsign <- ifelse(quaternion[[4]] > 0, 1, -1)
    nii$quatern_b <- quaternion[[1]] * qsign
    nii$quatern_c <- quaternion[[2]] * qsign
    nii$quatern_d <- quaternion[[3]] * qsign
    nii$qoffset_x <- vox2ras[1, 4]
    nii$qoffset_y <- vox2ras[2, 4]
    nii$qoffset_z <- vox2ras[3, 4]
    nii$srow_x <- vox2ras[1, ]
    nii$srow_y <- vox2ras[2, ]
    nii$srow_z <- vox2ras[3, ]


    RNifti::writeNifti(RNifti::asNifti(x, nii), file = con, ...)
  }

  get_xform_r <- function(con, useQuaternionFirst = TRUE) {
    img <- RNifti::readNifti(con, internal = TRUE, volumes = TRUE)
    xform <- RNifti::xform(img, useQuaternionFirst = useQuaternionFirst)
    xform
  }

  get_xform_py <- function(con) {
    stopifnot(file.exists(con))
    rpyants_module <- rpyANTs:::load_rpyants()
    img <- as_ANTsImage(con, strict = TRUE)
    xform <- rpyants_module$utils$transforms$get_xform(img)
    xform <- py_to_r(xform)
    xform
  }

  x <- array(rnorm(27), c(3,3,3))

  f <- tempfile(fileext = ".nii.gz")
  con <- f

  on.exit({ if(file.exists(f)) { unlink(f) } })


  # TODO: Make sure det(vox2ras) > 0
  vox2ras <- matrix(c(
    -0.345, 0, 0, 12,
    0, 0, -5.677, 15,
    0, -2.44, 0, 23,
    0, 0, 0, 1
  ), nrow = 4, ncol = 4, byrow = TRUE)
  if(det(vox2ras) < 0) {
    vox2ras[1, ] <- -vox2ras[1, ]
  }


  io_write_nii.array(x, f, vox2ras)
  xform_r <- get_xform_r(f)
  xform_py <- get_xform_py(f)
  testthat::expect_equal(solve(xform_r) %*% xform_py, diag(1, 4))


  vox2ras <- matrix(c(
    0, 0, -5.677, 15,
    -0.345, 0, 0, 12,
    0, -2.44, 0, 23,
    0, 0, 0, 1
  ), nrow = 4, ncol = 4, byrow = TRUE)
  if(det(vox2ras) < 0) {
    vox2ras[1, ] <- -vox2ras[1, ]
  }


  io_write_nii.array(x, f, vox2ras)
  xform_r <- get_xform_r(f)
  xform_py <- get_xform_py(f)
  testthat::expect_equal(solve(xform_r) %*% xform_py, diag(1, 4))


  vox2ras <- matrix(c(
    0, 0, -5.677, 15,
    -0.345, 0, 0, 12,
    0, -2.44, 0, 23,
    0, 0, 0, 1
  ), nrow = 4, ncol = 4, byrow = TRUE)
  vox2ras[1:3, ] <- vox2ras[1:3, ] + rnorm(12, sd = 0.01)
  if(det(vox2ras) < 0) {
    vox2ras[1, ] <- -vox2ras[1, ]
  }


  io_write_nii.array(x, f, vox2ras)
  xform_r <- get_xform_r(f)
  xform_py <- get_xform_py(f)
  testthat::expect_equal(solve(xform_r) %*% xform_py, diag(1, 4), tolerance = 1e-5)


})
