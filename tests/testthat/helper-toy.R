# Toy data for the tests: synthetic images that only take seconds to register.
# Nothing here reads real images.

# 2D image (96 x 96 mm) made of disks; each disk is c(row, column, radius,
# intensity) in mm. `spacing = 0.5` draws the same picture with twice the
# resolution
toy_disks <- function(..., spacing = 1) {
  n <- as.integer(round(96 / spacing))
  origin <- spacing / 2 - 0.5
  rows <- matrix(origin + spacing * (seq_len(n) - 1), n, n)
  cols <- t(rows)
  arr <- matrix(0, n, n)
  for (disk in list(...)) {
    arr <- arr + ((rows - disk[1])^2 + (cols - disk[2])^2 <= disk[3]^2) * disk[4]
  }
  img <- ants$from_numpy(
    np_array(arr, dtype = "float32"),
    origin = list(origin, origin), spacing = list(spacing, spacing)
  )
  ants$smooth_image(img, 1.5)
}

# 3D head made of ellipsoids (physical coordinates in mm, 2 mm voxels): scalp,
# skull, CSF, gray and white matter, ventricles, deep gray matter, brain stem,
# and cerebellum. Arguments `shift`, `scale`, `rotate` (degrees) and `ventricle`
# (relative size) make a different "subject"
toy_head <- function(shape, origin_offset = c(0, 0, 0), shift = c(0, 0, 0),
                     scale = 1, rotate = 0, ventricle = 1) {
  spacing <- c(2, 2, 2)
  origin <- -(shape - 1) * spacing / 2 + origin_offset
  grid <- expand.grid(
    x = origin[1] + spacing[1] * (seq_len(shape[1]) - 1) - shift[1],
    y = origin[2] + spacing[2] * (seq_len(shape[2]) - 1) - shift[2],
    z = origin[3] + spacing[3] * (seq_len(shape[3]) - 1) - shift[3]
  )
  theta <- rotate / 180 * pi
  x <- (cos(theta) * grid$x + sin(theta) * grid$y) / scale
  y <- (-sin(theta) * grid$x + cos(theta) * grid$y) / scale
  z <- grid$z / scale
  ellipsoid <- function(center, radius) {
    ((x - center[1]) / radius[1])^2 + ((y - center[2]) / radius[2])^2 +
      ((z - center[3]) / radius[3])^2 <= 1
  }

  brain <- ellipsoid(c(0, 0, 0), c(64, 78, 62))
  ventricles <- ellipsoid(c(-10, 0, 8), c(7, 22, 9) * ventricle) |
    ellipsoid(c(10, 0, 8), c(7, 22, 9) * ventricle)
  deep_gray <- ellipsoid(c(-24, 0, 4), c(9, 16, 9)) |
    ellipsoid(c(24, 0, 4), c(9, 16, 9))

  # Labels follow `deep_atropos`:
  #   1 CSF, 2 gray matter, 3 white matter, 4 deep gray matter, 5 brain stem,
  #   6 cerebellum
  seg <- integer(nrow(grid))
  seg[brain] <- 1L
  seg[ellipsoid(c(0, 0, 0), c(58, 72, 56))] <- 2L
  seg[ellipsoid(c(0, 0, 0), c(50, 64, 48))] <- 3L
  seg[brain & deep_gray] <- 4L
  seg[brain & ventricles] <- 1L
  seg[brain & ellipsoid(c(0, -10, -40), c(10, 10, 22))] <- 5L
  seg[brain & ellipsoid(c(0, -45, -38), c(34, 24, 18))] <- 6L

  t1 <- numeric(nrow(grid))
  t1[ellipsoid(c(0, 0, 0), c(78, 92, 78))] <- 60    # scalp
  t1[ellipsoid(c(0, 0, 0), c(70, 84, 70))] <- 12    # skull
  t1[seg > 0] <- c(30, 70, 100, 85, 92, 78)[seg[seg > 0]]

  as_image <- function(values, smooth = 0) {
    img <- ants$from_numpy(
      np_array(array(values, dim = shape), dtype = "float32"),
      origin = as.list(origin), spacing = as.list(spacing)
    )
    if (smooth > 0) {
      img <- ants$smooth_image(img, smooth)
    }
    img
  }

  list(t1 = t1, seg = seg, brain = brain, as_image = as_image)
}

# Writes a toy template folder (`T1.nii.gz`, `T1_brainmask.nii.gz`, and
# `atropos_{0..7}.nii.gz`: segmentation, background, then the six tissue
# probability maps) and a toy subject whose head is shifted, rotated, smaller,
# and has larger ventricles. The subject has a second image of the same head
# with another contrast (bright CSF, like a 'T2')
toy_head_dataset <- function(root) {
  template_dir <- file.path(root, "template")
  dir.create(template_dir, recursive = TRUE)

  template <- toy_head(c(88L, 104L, 88L))
  template$as_image(template$t1, smooth = 1.5)$to_file(
    file.path(template_dir, "T1.nii.gz"))
  template$as_image(template$brain)$to_file(
    file.path(template_dir, "T1_brainmask.nii.gz"))
  template$as_image(template$seg)$to_file(
    file.path(template_dir, "atropos_0.nii.gz"))
  for (label in 0:6) {
    template$as_image(template$seg == label, smooth = 1.5)$to_file(
      file.path(template_dir, sprintf("atropos_%d.nii.gz", label + 1L)))
  }

  subject <- toy_head(c(96L, 108L, 92L), origin_offset = c(3, -2, 1),
                      shift = c(4, -6, 3), scale = 0.94, rotate = 7,
                      ventricle = 1.3)
  # deterministic texture instead of random noise
  texture <- 1.5 * sin(seq_along(subject$t1) * 12.9898) * (subject$t1 > 0)
  subject$as_image(subject$t1 * 0.8 + texture, smooth = 1)$to_file(
    file.path(root, "subject_T1.nii.gz"))
  t2 <- c(0, 100, 60, 40, 50, 45, 55)[subject$seg + 1L] + (subject$t1 == 60) * 30
  subject$as_image(t2 + texture, smooth = 1)$to_file(
    file.path(root, "subject_T2.nii.gz"))
  subject$as_image(subject$seg)$to_file(file.path(root, "subject_seg.nii.gz"))

  list(
    template_t1 = file.path(template_dir, "T1.nii.gz"),
    template_seg = file.path(template_dir, "atropos_0.nii.gz"),
    subject_t1 = file.path(root, "subject_T1.nii.gz"),
    subject_t2 = file.path(root, "subject_T2.nii.gz"),
    subject_seg = file.path(root, "subject_seg.nii.gz")
  )
}

# Records every call to `ants.registration`: its `type_of_transform`, the
# weights of the additional metrics (`multivariate_extras`), whether a fixed
# mask is given, the settings of affine transforms if given (iterations, shrink
# factors, smoothing sigmas), the class of `initial_transform`, and whether
# `restrict_transformation` is given. The registration itself runs unchanged;
# call `restore()` when done
spy_ants_registration <- function() {
  env <- reticulate::py_run_string(paste(
    "def spy_on(registration, calls):",
    "    def join(values, sep):",
    "        return '' if values is None else sep.join(str(v) for v in values)",
    "    def spy(*args, **kwargs):",
    "        extras = kwargs.get('multivariate_extras', None)",
    "        calls.append({",
    "            'type': kwargs.get('type_of_transform', 'SyN'),",
    "            'weights': join(None if extras is None else [e[3] for e in extras], ', '),",
    "            'masked': kwargs.get('mask', None) is not None,",
    "            'affine': ' '.join(join(kwargs.get(key, None), 'x') for key in (",
    "                'aff_iterations', 'aff_shrink_factors', 'aff_smoothing_sigmas')).strip(),",
    "            'initial': type(kwargs.get('initial_transform', None)).__name__,",
    "            'restricted': kwargs.get('restrict_transformation', None) is not None,",
    "        })",
    "        return registration(*args, **kwargs)",
    "    return spy",
    sep = "\n"
  ), local = TRUE, convert = FALSE)

  ants_module <- import("ants", convert = FALSE)
  registration <- ants_module$registration
  calls <- r_to_py(list())
  py_set_attr(ants_module, "registration", env$spy_on(registration, calls))

  list(
    # data frame with one row per call, in the order of the calls
    calls = function() {
      do.call(rbind, lapply(py_to_r(calls), as.data.frame))
    },
    restore = function() {
      py_set_attr(ants_module, "registration", registration)
    }
  )
}

# Makes 'Python', hence 'ANTsPy', create its temporary files in `dir`, so that
# a test can see and remove them. Returns a function that restores the
# previous folder
use_python_tempdir <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tempfile_module <- import("tempfile", convert = FALSE)
  previous <- tempfile_module$tempdir
  py_set_attr(tempfile_module, "tempdir", dir)
  function() {
    py_set_attr(tempfile_module, "tempdir", previous)
  }
}

# Evaluates `expr` while `os.listdir` of 'Python' returns sorted names; the
# order is otherwise up to the file system
with_sorted_listdir <- function(expr) {
  env <- reticulate::py_run_string(paste(
    "def sort_listing(listdir):",
    "    return lambda *args, **kwargs: sorted(listdir(*args, **kwargs))",
    sep = "\n"
  ), local = TRUE, convert = FALSE)

  os <- import("os", convert = FALSE)
  listdir <- os$listdir
  py_set_attr(os, "listdir", env$sort_listing(listdir))
  on.exit({ py_set_attr(os, "listdir", listdir) })
  expr
}

# Overlap (Dice) of each tissue label between two segmentation arrays
toy_dice <- function(seg1, seg2, labels = 1:6) {
  vapply(labels, function(label) {
    2 * sum(seg1 == label & seg2 == label) /
      (sum(seg1 == label) + sum(seg2 == label))
  }, 0.0)
}
