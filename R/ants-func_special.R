# Special treatments

#' @rdname ants.invariant_image_similarity
#'
#' @export
ants.invariant_image_similarity <- function(
    image1, image2, local_search_iterations = 0L, metric = "MI",
    thetas = c( 0.,  90., 180., 270., 360. ),
    thetas2 = c( 0.,  90., 180., 270., 360. ),
    thetas3 = c( 0.,  90., 180., 270., 360. ),
    scale_image = 1L, do_reflection = FALSE, txfn = NULL, transform = "Affine") {
  ants <- load_ants()

  thetas <- np_array(thetas)
  thetas2 <- np_array(thetas2)
  thetas3 <- np_array(thetas3)

  ants$invariant_image_similarity(
    image1 = image1,
    image2 = image2,
    local_search_iterations = local_search_iterations,
    metric = metric,
    thetas = thetas,
    thetas2 = thetas2,
    thetas3 = thetas3,
    scale_image = scale_image,
    do_reflection = do_reflection,
    txfn = txfn,
    transform = transform
  )
}



#' @rdname ants.labels_to_matrix
#'
#' @export
ants.labels_to_matrix <- function(image, mask, target_labels = NULL, missing_val = NaN) {
  ants <- load_ants()
  ants$labels_to_matrix(
    image = image,
    mask = mask,
    target_labels = target_labels,
    missing_val = missing_val
  )
}


#' @rdname ants.sparse_decom2
#'
#' @export
ants.sparse_decom2 <- function(inmatrix, inmask = tuple(NULL, NULL), sparseness = tuple(0.01, 0.01), nvecs = 3L, its = 20L, cthresh = tuple(0L, 0L), statdir = NULL, perms = 0L, uselong = 0L, z = 0L, smooth = 0L, robust = 0L, mycoption = 0L, initialization_list = list(), initialization_list2 = list(), ell1 = 10L, prior_weight = 0L, verbose = FALSE, rejector = 0L, max_based = FALSE, version = 1L) {
  ants <- load_ants()
  ants$sparse_decom2(
    inmatrix = inmatrix,
    inmask = inmask,
    sparseness = sparseness,
    nvecs = nvecs,
    its = its,
    cthresh = cthresh,
    statdir = statdir,
    perms = perms,
    uselong = uselong,
    z = z,
    smooth = smooth,
    robust = robust,
    mycoption = mycoption,
    initialization_list = initialization_list,
    initialization_list2 = initialization_list2,
    ell1 = ell1,
    prior_weight = prior_weight,
    verbose = verbose,
    rejector = rejector,
    max_based = max_based,
    version = version
  )
}


#' @rdname ants.apply_ants_transform
#'
#' @export
ants.apply_ants_transform <- function(transform, data, data_type = "point", reference = NULL, ...) {
  ants <- load_ants()
  ants$apply_ants_transform(
    transform = transform,
    data = data,
    data_type = data_type,
    reference = reference,
    ...
  )
}
