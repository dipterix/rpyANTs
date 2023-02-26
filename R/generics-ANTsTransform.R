# DIPSAUS DEBUG START
# x <- ants$new_ants_transform()
# fi <- ants$image_read(ants$get_ants_data('r16'))
# mo <- ants$image_read(ants$get_ants_data('r64'))
#
# # resample to speed up this example
# fi <- ants$resample_image(fi, list(60L,60L), TRUE, 0L)
# mo <- ants$resample_image(mo, list(60L,60L), TRUE, 0L)
#
# # SDR transform
# transform <- ants_registration(
#   fixed=fi, moving=mo, type_of_transform = 'SyN' )
# y <- ants$image_read(transform$fwdtransforms[0])
# x <- ants$transform_from_displacement_field(y)

#' @title Check if an object is a 3D 'affine' transform matrix
#' @param x R or Python object, accepted forms are numeric \code{matrix},
#' \code{'ANTsTransform'}, or \code{character} (path to transform matrix)
#' @param strict whether the last element should be always 1
#' @param ... passed to other methods
#' @returns A logical value whether the object can be loaded as a 4-by-4 matrix.
#'
#' @examples
#'
#' # not affine
#' is_affine3D(1)
#'
#' # 3x3 matrix is not as it is treated as 2D transform
#' is_affine3D(matrix(rnorm(9), nrow = 3))
#'
#' # 3x4 matrix
#' x <- matrix(rnorm(12), nrow = 3)
#' is_affine3D(x)
#'
#' # 4x4 matrix
#' x <- rbind(x, c(0,0,0,1))
#' is_affine3D(x)
#'
#' if(interactive() && ants_available()) {
#'
#'   ants <- load_ants()
#'   x <- ants$new_ants_transform(dimension = 3L)
#'   is_affine3D(x)
#'
#'   # save the parameters
#'   f <- tempfile(fileext = ".mat")
#'   ants$write_transform(x, f)
#'   is_affine3D(f)
#'
#' }
#'
#'
#'
#'
#' @export
is_affine3D <- function(x, ...) {
  UseMethod("is_affine3D")
}

#' @rdname is_affine3D
#' @export
is_affine3D.default <- function(x, strict = TRUE, ...) {
  if(!is.matrix(x)) {
    re <- FALSE
    if(is.character(x)) {
      re <- tryCatch({
        ants <- load_ants()
        x <- ants$read_transform(x)
        is_affine3D.ants.core.ants_transform.ANTsTransform(x)
      }, error = function(e){
        FALSE
      })
    }
    return(re)
  }
  if(!is.numeric(x)) { return(FALSE) }
  if(ncol(x) != 4) { return(FALSE) }
  if(nrow(x) == 4) {
    if(strict && x[[16]] != 1) {
      return(FALSE)
    }
    return(TRUE)
  }
  if(nrow(x) == 3) { return(TRUE) }
  return(TRUE)
}

#' @rdname is_affine3D
#' @export
is_affine3D.ants.core.ants_transform.ANTsTransform <- function(x, ...) {
  if(!isTRUE(to_r(x$dimension) == 3)) {
    return(FALSE)
  }
  if(is_py_inherits(x$transform_type)) {
    transform_type <- to_r(x$transform_type)
  } else {
    transform_type <- x$transform_type
  }

  if(!transform_type %in% AFFINE_TRANSFORM_TYPES) { return(FALSE) }
  return(TRUE)
}

#' @title Convert to \code{'ANTsTransform'}
#' @param x 'affine' matrix or \code{'numpy'} array, character path to the
#' matrix, \code{'ANTsTransform'}, \code{'ANTsImage'} as displacement field.
#' @param dimension expected transform space dimension; default is 3
#' @param ... passed to other methods
#' @returns An \code{'ANTsTransform'} object
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   mat <- matrix(c(
#'     0, -1, 0, 128,
#'     1, 0, 0, -128,
#'     0, 0, -1, 128,
#'     0, 0,  0,   1
#'   ), ncol = 4, byrow = TRUE)
#'
#'   trans <- as_ANTsTransform(mat)
#'   trans
#'
#'   # apply transform
#'   trans$apply_to_point(c(120, 400, 1))
#'
#'   # same results
#'   mat %*% c(120, 400, 1, 1)
#'
#'   trans[] == mat
#'
#' }
#'
#' @export
as_ANTsTransform <- function(x, ...) {
  UseMethod("as_ANTsTransform")
}

#' @rdname as_ANTsTransform
#' @export
as_ANTsTransform.default <- function(x, dimension = 3, ...) {
  dm <- dim(x)

  is_affine <- TRUE
  if(is_py_inherits(dimension)) {
    dimension <- to_r(dimension)
  }
  dimension <- as.integer(dimension)[[1]]

  if(length(dm) != 2) {
    is_affine <- FALSE
  } else if(dm[[1]] > dimension+1 || dm[[2]] > dimension+1) {
    is_affine <- FALSE
  }

  if(!is_affine) {
    stop("as_ANTsTransform: non-linear transform is not yet supported.")
    # if(dm[[1]] != dimension) {
    #   stop(sprintf("as_ANTsTransform: is this a %sD non-linear transform? Please specify `dimension` explicitly", dm[[1]]))
    # }
    # if(length(dm) != dimension + 1) {
    #   stop("as_ANTsTransform: input array does not seem to be a displacement field nor affine transform.")
    # }
  }

  ants <- load_ants()

  if( is_affine ) {
    y <- cbind(diag(rep(1, dimension)), 0)
    dm[[1]] <- min(dm[[1]], dimension)
    dm[[2]] <- min(dm[[2]], dimension + 1)
    y[seq_len(dm[[1]]), seq_len(dm[[2]])] <- x[seq_len(dm[[1]]), seq_len(dm[[2]])]

    return(
      ants$create_ants_transform(
        transform_type = "AffineTransform",
        dimension = dimension,
        matrix = y,
        translation = y[, dimension + 1]
      )
    )
  }
#
#   ants$create_ants_transform(
#     transform_type = "DisplacementFieldTransform",
#     dimension = dimension, parameters = as.vector(x)
#   )
#   # non-linear
#   np <- import("numpy", convert = FALSE)
#   z <- np$asarray(x)
#   ants$ANTsImage(dimension = dimension, components = )
#   field <- ants$from_numpy(data = z, direction = diag(rep(1.0, dimension)), has_components = TRUE, is_rgb = FALSE)
#   return(ants$transform_from_displacement_field(field = field))
}

#' @rdname as_ANTsTransform
#' @export
as_ANTsTransform.ants.core.ants_transform.ANTsTransform <- function(x, ...) {
  return(x)
}

#' @rdname as_ANTsTransform
#' @export
as_ANTsTransform.ants.core.ants_image.ANTsImage <- function(x, ...) {
  ants <- load_ants()
  return(ants$transform_from_displacement_field(field = x))
}

#' @rdname as_ANTsTransform
#' @export
as_ANTsTransform.numpy.ndarray <- function(x, ...) {
  as_ANTsTransform.default(to_r(x), ...)
}

#' @rdname as_ANTsTransform
#' @export
as_ANTsTransform.character <- function(x, ...) {
  ants <- load_ants()
  ants$read_transform(filename = x, ...)
}

#' @export
as.matrix.ants.core.ants_transform.ANTsTransform <- function(x, ...) {

  if(is_py_inherits(x$transform_type)) {
    transform_type <- to_r(x$transform_type)
  } else {
    transform_type <- x$transform_type
  }

  if(!isTRUE(transform_type %in% AFFINE_TRANSFORM_TYPES)) {
    stop("This ANTsTransform is not affine/linear. Cannot convert to matrix")
  }
  ndims <- to_r(x$dimension)

  pfrom <- diag(rep(1.0, ndims))
  tx <- apply(pfrom, 2, function(v) {
    c(as.double(to_r(x$apply_to_vector(v))), 0)
  })

  trans <- c(as.double(to_r(x$apply_to_point(rep(0.0, ndims)))), 1)
  unname(cbind(tx, trans))
}

#' @export
as.array.ants.core.ants_transform.ANTsTransform <- function(x, displacement_field = NULL, ...) {
  if(is_py_inherits(x$transform_type)) {
    transform_type <- to_r(x$transform_type)
  } else {
    transform_type <- x$transform_type
  }


  if(transform_type %in% AFFINE_TRANSFORM_TYPES) {
    return(as.matrix(x, ...))
  }

  if(!transform_type %in% "DisplacementFieldTransform") {
    stop("as.array: Unsupported ANTsTransform type: ", transform_type)
  }

  if(!is_py_inherits(displacement_field, "ANTsImage")) {
    stop("as.array: `displacement_field` must be an ANTsImage")
  }

  ants <- load_ants()
  b <- ants$transform_to_displacement_field(x, displacement_field$clone())
  return(as.array(b))
}

#' @export
`[.ants.core.ants_transform.ANTsTransform` <- function(
    x, ..., drop = TRUE, displacement_field = NULL) {
  as.array(x, displacement_field = displacement_field)[..., drop = drop]
}

#' @export
dim.ants.core.ants_transform.ANTsTransform <- function(x) {
  if(is_py_inherits(x$transform_type)) {
    transform_type <- to_r(x$transform_type)
  } else {
    transform_type <- x$transform_type
  }

  if(!isTRUE(transform_type %in% AFFINE_TRANSFORM_TYPES)) {
    stop("Cannot obtain dim(x) for non-linear transform")
  }
  ndims <- to_r(x$dimension) + 1
  return(c(ndims, ndims))
}
