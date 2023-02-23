

#' @rdname ants.abp_n4
#'
#' @export
ants.abp_n4 <- function(image, intensity_truncation = list(0.025, 0.975, 256L), mask = NULL, usen3 = FALSE) {
  ants <- load_ants()
  ants$abp_n4(
    image = image,
    intensity_truncation = intensity_truncation,
    mask = mask,
    usen3 = usen3
  )
}


#' @rdname ants.add_noise_to_image
#'
#' @export
ants.add_noise_to_image <- function(image, noise_model, noise_parameters) {
  ants <- load_ants()
  ants$add_noise_to_image(
    image = image,
    noise_model = noise_model,
    noise_parameters = noise_parameters
  )
}


#' @rdname ants.affine_initializer
#'
#' @export
ants.affine_initializer <- function(fixed_image, moving_image, search_factor = 20L, radian_fraction = 0.1, use_principal_axis = FALSE, local_search_iterations = 10L, mask = NULL, txfn = NULL) {
  ants <- load_ants()
  ants$affine_initializer(
    fixed_image = fixed_image,
    moving_image = moving_image,
    search_factor = search_factor,
    radian_fraction = radian_fraction,
    use_principal_axis = use_principal_axis,
    local_search_iterations = local_search_iterations,
    mask = mask,
    txfn = txfn
  )
}


#' @rdname ants.allclose
#'
#' @export
ants.allclose <- function(image1, image2) {
  ants <- load_ants()
  ants$allclose(
    image1 = image1,
    image2 = image2
  )
}


#' @rdname ants.anti_alias
#'
#' @export
ants.anti_alias <- function(image) {
  ants <- load_ants()
  ants$anti_alias(
    image = image
  )
}


#' @rdname ants.ANTsImage
#'
#' @export
ants.ANTsImage <- function(pixeltype = "float", dimension = 3L, components = 1L, pointer = NULL, is_rgb = FALSE, label_image = NULL) {
  ants <- load_ants()
  ants$ANTsImage(
    pixeltype = pixeltype,
    dimension = dimension,
    components = components,
    pointer = pointer,
    is_rgb = is_rgb,
    label_image = label_image
  )
}


#' @rdname ants.ANTsTransform
#'
#' @export
ants.ANTsTransform <- function(precision = "float", dimension = 3L, transform_type = "AffineTransform", pointer = NULL) {
  ants <- load_ants()
  ants$ANTsTransform(
    precision = precision,
    dimension = dimension,
    transform_type = transform_type,
    pointer = pointer
  )
}





#' @rdname ants.apply_ants_transform
#'
#' @export
ants.apply_ants_transform_to_image <- function(transform, image, reference, interpolation = "linear") {
  ants <- load_ants()
  ants$apply_ants_transform_to_image(
    transform = transform,
    image = image,
    reference = reference,
    interpolation = interpolation
  )
}


#' @rdname ants.apply_ants_transform
#'
#' @export
ants.apply_ants_transform_to_point <- function(transform, point) {
  ants <- load_ants()
  ants$apply_ants_transform_to_point(
    transform = transform,
    point = point
  )
}


#' @rdname ants.apply_ants_transform
#'
#' @export
ants.apply_ants_transform_to_vector <- function(transform, vector) {
  ants <- load_ants()
  ants$apply_ants_transform_to_vector(
    transform = transform,
    vector = vector
  )
}


#' @rdname ants.apply_transforms
#'
#' @export
ants.apply_transforms <- function(fixed, moving, transformlist, interpolator = "linear", imagetype = 0L, whichtoinvert = NULL, compose = NULL, defaultvalue = 0L, verbose = FALSE) {
  ants <- load_ants()
  ants$apply_transforms(
    fixed = fixed,
    moving = moving,
    transformlist = transformlist,
    interpolator = interpolator,
    imagetype = imagetype,
    whichtoinvert = whichtoinvert,
    compose = compose,
    defaultvalue = defaultvalue,
    verbose = verbose
  )
}


#' @rdname ants.apply_ants_transform
#'
#' @export
ants.apply_transforms_to_points <- function(dim, points, transformlist, whichtoinvert = NULL, verbose = FALSE) {
  ants <- load_ants()
  ants$apply_transforms_to_points(
    dim = dim,
    points = points,
    transformlist = transformlist,
    whichtoinvert = whichtoinvert,
    verbose = verbose
  )
}


#' @rdname ants.atropos
#'
#' @export
ants.atropos <- function(a, x, i = "Kmeans[3]", m = "[0.2,1x1]", c = "[5,0]", priorweight = 0.25) {
  ants <- load_ants()
  ants$atropos(
    a = a,
    x = x,
    i = i,
    m = m,
    c = c,
    priorweight = priorweight
  )
}


#' @rdname ants.bandpass_filter_matrix
#'
#' @export
ants.bandpass_filter_matrix <- function(matrix, tr = 1L, lowf = 0.01, highf = 0.1, order = 3L) {
  ants <- load_ants()
  ants$bandpass_filter_matrix(
    matrix = matrix,
    tr = tr,
    lowf = lowf,
    highf = highf,
    order = order
  )
}


#' @rdname ants.build_template
#'
#' @export
ants.build_template <- function(initial_template = NULL, image_list = NULL, iterations = 3L, gradient_step = 0.2, blending_weight = 0.75, weights = NULL) {
  ants <- load_ants()
  ants$build_template(
    initial_template = initial_template,
    image_list = image_list,
    iterations = iterations,
    gradient_step = gradient_step,
    blending_weight = blending_weight,
    weights = weights
  )
}


#' @rdname ants.compcor
#'
#' @export
ants.compcor <- function(boldImage, ncompcor = 4L, quantile = 0.975, mask = NULL, filter_type = FALSE, degree = 2L) {
  ants <- load_ants()
  ants$compcor(
    boldImage = boldImage,
    ncompcor = ncompcor,
    quantile = quantile,
    mask = mask,
    filter_type = filter_type,
    degree = degree
  )
}


#' @rdname ants.compose_ants_transforms
#'
#' @export
ants.compose_ants_transforms <- function(transform_list) {
  ants <- load_ants()
  ants$compose_ants_transforms(
    transform_list = transform_list
  )
}


#' @rdname ants.compose_displacement_fields
#'
#' @export
ants.compose_displacement_fields <- function(displacement_field, warping_field) {
  ants <- load_ants()
  ants$compose_displacement_fields(
    displacement_field = displacement_field,
    warping_field = warping_field
  )
}


#' @rdname ants.convolve_image
#'
#' @export
ants.convolve_image <- function(image, kernel_image, crop = TRUE) {
  ants <- load_ants()
  ants$convolve_image(
    image = image,
    kernel_image = kernel_image,
    crop = crop
  )
}


#' @rdname ants.copy_image_info
#'
#' @export
ants.copy_image_info <- function(reference, target) {
  ants <- load_ants()
  ants$copy_image_info(
    reference = reference,
    target = target
  )
}


#' @rdname ants.create_ants_metric
#'
#' @export
ants.create_ants_metric <- function(fixed, moving, metric_type = "MeanSquares", fixed_mask = NULL, moving_mask = NULL, sampling_strategy = "regular", sampling_percentage = 1L) {
  ants <- load_ants()
  ants$create_ants_metric(
    fixed = fixed,
    moving = moving,
    metric_type = metric_type,
    fixed_mask = fixed_mask,
    moving_mask = moving_mask,
    sampling_strategy = sampling_strategy,
    sampling_percentage = sampling_percentage
  )
}


#' @rdname ants.create_ants_transform
#'
#' @export
ants.create_ants_transform <- function(transform_type = "AffineTransform", precision = "float", dimension = 3L, matrix = NULL, offset = NULL, center = NULL, translation = NULL, parameters = NULL, fixed_parameters = NULL, displacement_field = NULL, supported_types = FALSE) {
  ants <- load_ants()
  ants$create_ants_transform(
    transform_type = transform_type,
    precision = precision,
    dimension = dimension,
    matrix = matrix,
    offset = offset,
    center = center,
    translation = translation,
    parameters = parameters,
    fixed_parameters = fixed_parameters,
    displacement_field = displacement_field,
    supported_types = supported_types
  )
}


#' @rdname ants.create_jacobian_determinant_image
#'
#' @export
ants.create_jacobian_determinant_image <- function(domain_image, tx, do_log = FALSE, geom = FALSE) {
  ants <- load_ants()
  ants$create_jacobian_determinant_image(
    domain_image = domain_image,
    tx = tx,
    do_log = do_log,
    geom = geom
  )
}


#' @rdname ants.create_tiled_mosaic
#'
#' @export
ants.create_tiled_mosaic <- function(image, rgb = NULL, mask = NULL, overlay = NULL, output = NULL, alpha = 1.0, direction = 0L, pad_or_crop = NULL, slices = NULL, flip_slice = NULL, permute_axes = FALSE) {
  ants <- load_ants()
  ants$create_tiled_mosaic(
    image = image,
    rgb = rgb,
    mask = mask,
    overlay = overlay,
    output = output,
    alpha = alpha,
    direction = direction,
    pad_or_crop = pad_or_crop,
    slices = slices,
    flip_slice = flip_slice,
    permute_axes = permute_axes
  )
}


#' @rdname ants.create_warped_grid
#'
#' @export
ants.create_warped_grid <- function(image, grid_step = 10L, grid_width = 2L, grid_directions = list(TRUE, TRUE), fixed_reference_image = NULL, transform = NULL, foreground = 1L, background = 0L) {
  ants <- load_ants()
  ants$create_warped_grid(
    image = image,
    grid_step = grid_step,
    grid_width = grid_width,
    grid_directions = grid_directions,
    fixed_reference_image = fixed_reference_image,
    transform = transform,
    foreground = foreground,
    background = background
  )
}


#' @rdname ants.crop_image
#'
#' @export
ants.crop_image <- function(image, label_image = NULL, label = 1L) {
  ants <- load_ants()
  ants$crop_image(
    image = image,
    label_image = label_image,
    label = label
  )
}


#' @rdname ants.crop_indices
#'
#' @export
ants.crop_indices <- function(image, lowerind, upperind) {
  ants <- load_ants()
  ants$crop_indices(
    image = image,
    lowerind = lowerind,
    upperind = upperind
  )
}


#' @rdname ants.decrop_image
#'
#' @export
ants.decrop_image <- function(cropped_image, full_image) {
  ants <- load_ants()
  ants$decrop_image(
    cropped_image = cropped_image,
    full_image = full_image
  )
}


#' @rdname ants.denoise_image
#'
#' @export
ants.denoise_image <- function(image, mask = NULL, shrink_factor = 1L, p = 1L, r = 3L, noise_model = "Rician", v = 0L) {
  ants <- load_ants()
  ants$denoise_image(
    image = image,
    mask = mask,
    shrink_factor = shrink_factor,
    p = p,
    r = r,
    noise_model = noise_model,
    v = v
  )
}


#' @rdname ants.dicom_read
#'
#' @export
ants.dicom_read <- function(directory, pixeltype = "float") {
  ants <- load_ants()
  ants$dicom_read(
    directory = directory,
    pixeltype = pixeltype
  )
}


#' @rdname ants.eig_seg
#'
#' @export
ants.eig_seg <- function(mask, img_list, apply_segmentation_to_images = FALSE, cthresh = 0L, smooth = 1L) {
  ants <- load_ants()
  ants$eig_seg(
    mask = mask,
    img_list = img_list,
    apply_segmentation_to_images = apply_segmentation_to_images,
    cthresh = cthresh,
    smooth = smooth
  )
}


#' @rdname ants.fit_bspline_displacement_field
#'
#' @export
ants.fit_bspline_displacement_field <- function(displacement_field = NULL, displacement_weight_image = NULL, displacement_origins = NULL, displacements = NULL, displacement_weights = NULL, origin = NULL, spacing = NULL, size = NULL, direction = NULL, number_of_fitting_levels = 4L, mesh_size = 1L, spline_order = 3L, enforce_stationary_boundary = TRUE, estimate_inverse = FALSE, rasterize_points = FALSE) {
  ants <- load_ants()
  ants$fit_bspline_displacement_field(
    displacement_field = displacement_field,
    displacement_weight_image = displacement_weight_image,
    displacement_origins = displacement_origins,
    displacements = displacements,
    displacement_weights = displacement_weights,
    origin = origin,
    spacing = spacing,
    size = size,
    direction = direction,
    number_of_fitting_levels = number_of_fitting_levels,
    mesh_size = mesh_size,
    spline_order = spline_order,
    enforce_stationary_boundary = enforce_stationary_boundary,
    estimate_inverse = estimate_inverse,
    rasterize_points = rasterize_points
  )
}


#' @rdname ants.fit_bspline_object_to_scattered_data
#'
#' @export
ants.fit_bspline_object_to_scattered_data <- function(scattered_data, parametric_data, parametric_domain_origin, parametric_domain_spacing, parametric_domain_size, is_parametric_dimension_closed = NULL, data_weights = NULL, number_of_fitting_levels = 4L, mesh_size = 1L, spline_order = 3L) {
  ants <- load_ants()
  ants$fit_bspline_object_to_scattered_data(
    scattered_data = scattered_data,
    parametric_data = parametric_data,
    parametric_domain_origin = parametric_domain_origin,
    parametric_domain_spacing = parametric_domain_spacing,
    parametric_domain_size = parametric_domain_size,
    is_parametric_dimension_closed = is_parametric_dimension_closed,
    data_weights = data_weights,
    number_of_fitting_levels = number_of_fitting_levels,
    mesh_size = mesh_size,
    spline_order = spline_order
  )
}


#' @rdname ants.fit_transform_to_paired_points
#'
#' @export
ants.fit_transform_to_paired_points <- function(moving_points, fixed_points, transform_type = "affine", regularization = 1e-06, domain_image = NULL, number_of_fitting_levels = 4L, mesh_size = 1L, spline_order = 3L, enforce_stationary_boundary = TRUE, displacement_weights = NULL, number_of_compositions = 10L, composition_step_size = 0.5, sigma = 0.0, convergence_threshold = 0.0, number_of_integration_points = 2L, rasterize_points = FALSE, verbose = FALSE) {
  ants <- load_ants()
  ants$fit_transform_to_paired_points(
    moving_points = moving_points,
    fixed_points = fixed_points,
    transform_type = transform_type,
    regularization = regularization,
    domain_image = domain_image,
    number_of_fitting_levels = number_of_fitting_levels,
    mesh_size = mesh_size,
    spline_order = spline_order,
    enforce_stationary_boundary = enforce_stationary_boundary,
    displacement_weights = displacement_weights,
    number_of_compositions = number_of_compositions,
    composition_step_size = composition_step_size,
    sigma = sigma,
    convergence_threshold = convergence_threshold,
    number_of_integration_points = number_of_integration_points,
    rasterize_points = rasterize_points,
    verbose = verbose
  )
}


#' @rdname ants.from_nibabel
#'
#' @export
ants.from_nibabel <- function(nib_image) {
  ants <- load_ants()
  ants$from_nibabel(
    nib_image = nib_image
  )
}


#' @rdname ants.from_numpy
#'
#' @export
ants.from_numpy <- function(data, origin = NULL, spacing = NULL, direction = NULL, has_components = FALSE, is_rgb = FALSE) {
  ants <- load_ants()
  ants$from_numpy(
    data = data,
    origin = origin,
    spacing = spacing,
    direction = direction,
    has_components = has_components,
    is_rgb = is_rgb
  )
}


#' @rdname ants.fsl2antstransform
#'
#' @export
ants.fsl2antstransform <- function(matrix, reference, moving) {
  ants <- load_ants()
  ants$fsl2antstransform(
    matrix = matrix,
    reference = reference,
    moving = moving
  )
}


#' @rdname ants.functional_lung_segmentation
#'
#' @export
ants.functional_lung_segmentation <- function(image, mask = NULL, number_of_iterations = 2L, number_of_atropos_iterations = 5L, mrf_parameters = "[0.7,2x2x2]", number_of_clusters = 6L, cluster_centers = NULL, bias_correction = "n4", verbose = TRUE) {
  ants <- load_ants()
  ants$functional_lung_segmentation(
    image = image,
    mask = mask,
    number_of_iterations = number_of_iterations,
    number_of_atropos_iterations = number_of_atropos_iterations,
    mrf_parameters = mrf_parameters,
    number_of_clusters = number_of_clusters,
    cluster_centers = cluster_centers,
    bias_correction = bias_correction,
    verbose = verbose
  )
}


#' @rdname ants.fuzzy_spatial_cmeans_segmentation
#'
#' @export
ants.fuzzy_spatial_cmeans_segmentation <- function(image, mask = NULL, number_of_clusters = 4L, m = 2L, p = 1L, q = 1L, radius = 2L, max_number_of_iterations = 20L, convergence_threshold = 0.02, verbose = FALSE) {
  ants <- load_ants()
  ants$fuzzy_spatial_cmeans_segmentation(
    image = image,
    mask = mask,
    number_of_clusters = number_of_clusters,
    m = m,
    p = p,
    q = q,
    radius = radius,
    max_number_of_iterations = max_number_of_iterations,
    convergence_threshold = convergence_threshold,
    verbose = verbose
  )
}


#' @rdname ants.get_ants_data
#'
#' @export
ants.get_ants_data <- function(file_id = NULL, target_file_name = NULL, antsx_cache_directory = NULL) {
  ants <- load_ants()
  ants$get_ants_data(
    file_id = file_id,
    target_file_name = target_file_name,
    antsx_cache_directory = antsx_cache_directory
  )
}


#' @rdname ants.get_ants_transform_fixed_parameters
#'
#' @export
ants.get_ants_transform_fixed_parameters <- function(transform) {
  ants <- load_ants()
  ants$get_ants_transform_fixed_parameters(
    transform = transform
  )
}


#' @rdname ants.get_ants_transform_parameters
#'
#' @export
ants.get_ants_transform_parameters <- function(transform) {
  ants <- load_ants()
  ants$get_ants_transform_parameters(
    transform = transform
  )
}


#' @rdname ants.get_average_of_timeseries
#'
#' @export
ants.get_average_of_timeseries <- function(image, idx = NULL) {
  ants <- load_ants()
  ants$get_average_of_timeseries(
    image = image,
    idx = idx
  )
}


#' @rdname ants.get_canonical_views
#'
#' @export
ants.get_canonical_views <- function() {
  ants <- load_ants()
  ants$get_canonical_views(
)
}


#' @rdname ants.get_center_of_mass
#'
#' @export
ants.get_center_of_mass <- function(image) {
  ants <- load_ants()
  ants$get_center_of_mass(
    image = image
  )
}


#' @rdname ants.get_centroids
#'
#' @export
ants.get_centroids <- function(image, clustparam = 0L) {
  ants <- load_ants()
  ants$get_centroids(
    image = image,
    clustparam = clustparam
  )
}


#' @rdname ants.get_data
#'
#' @export
ants.get_data <- function(file_id = NULL, target_file_name = NULL, antsx_cache_directory = NULL) {
  ants <- load_ants()
  ants$get_data(
    file_id = file_id,
    target_file_name = target_file_name,
    antsx_cache_directory = antsx_cache_directory
  )
}


#' @rdname ants.get_direction
#'
#' @export
ants.get_direction <- function(image) {
  ants <- load_ants()
  ants$get_direction(
    image = image
  )
}


#' @rdname ants.get_lib_fn
#'
#' @export
ants.get_lib_fn <- function(string) {
  ants <- load_ants()
  ants$get_lib_fn(
    string = string
  )
}


#' @rdname ants.get_mask
#'
#' @export
ants.get_mask <- function(image, low_thresh = NULL, high_thresh = NULL, cleanup = 2L) {
  ants <- load_ants()
  ants$get_mask(
    image = image,
    low_thresh = low_thresh,
    high_thresh = high_thresh,
    cleanup = cleanup
  )
}


#' @rdname ants.get_neighborhood_at_voxel
#'
#' @export
ants.get_neighborhood_at_voxel <- function(image, center, kernel, physical_coordinates = FALSE) {
  ants <- load_ants()
  ants$get_neighborhood_at_voxel(
    image = image,
    center = center,
    kernel = kernel,
    physical_coordinates = physical_coordinates
  )
}


#' @rdname ants.get_neighborhood_in_mask
#'
#' @export
ants.get_neighborhood_in_mask <- function(image, mask, radius, physical_coordinates = FALSE, boundary_condition = NULL, spatial_info = FALSE, get_gradient = FALSE) {
  ants <- load_ants()
  ants$get_neighborhood_in_mask(
    image = image,
    mask = mask,
    radius = radius,
    physical_coordinates = physical_coordinates,
    boundary_condition = boundary_condition,
    spatial_info = spatial_info,
    get_gradient = get_gradient
  )
}


#' @rdname ants.get_orientation
#'
#' @export
ants.get_orientation <- function(image) {
  ants <- load_ants()
  ants$get_orientation(
    image = image
  )
}


#' @rdname ants.get_origin
#'
#' @export
ants.get_origin <- function(image) {
  ants <- load_ants()
  ants$get_origin(
    image = image
  )
}


#' @rdname ants.get_pointer_string
#'
#' @export
ants.get_pointer_string <- function(image) {
  ants <- load_ants()
  ants$get_pointer_string(
    image = image
  )
}


#' @rdname ants.get_possible_orientations
#'
#' @export
ants.get_possible_orientations <- function() {
  ants <- load_ants()
  ants$get_possible_orientations(
)
}


#' @rdname ants.get_spacing
#'
#' @export
ants.get_spacing <- function(image) {
  ants <- load_ants()
  ants$get_spacing(
    image = image
  )
}


#' @rdname ants.hausdorff_distance
#'
#' @export
ants.hausdorff_distance <- function(image1, image2) {
  ants <- load_ants()
  ants$hausdorff_distance(
    image1 = image1,
    image2 = image2
  )
}


#' @rdname ants.histogram_match_image
#'
#' @export
ants.histogram_match_image <- function(source_image, reference_image, number_of_histogram_bins = 255L, number_of_match_points = 64L, use_threshold_at_mean_intensity = FALSE) {
  ants <- load_ants()
  ants$histogram_match_image(
    source_image = source_image,
    reference_image = reference_image,
    number_of_histogram_bins = number_of_histogram_bins,
    number_of_match_points = number_of_match_points,
    use_threshold_at_mean_intensity = use_threshold_at_mean_intensity
  )
}


#' @rdname ants.ilr
#'
#' @export
ants.ilr <- function(data_frame, voxmats, ilr_formula, verbose = FALSE) {
  ants <- load_ants()
  ants$ilr(
    data_frame = data_frame,
    voxmats = voxmats,
    ilr_formula = ilr_formula,
    verbose = verbose
  )
}


#' @rdname ants.image_clone
#'
#' @export
ants.image_clone <- function(image, pixeltype = NULL) {
  ants <- load_ants()
  ants$image_clone(
    image = image,
    pixeltype = pixeltype
  )
}


#' @rdname ants.image_header_info
#'
#' @export
ants.image_header_info <- function(filename) {
  ants <- load_ants()
  ants$image_header_info(
    filename = filename
  )
}


#' @rdname ants.image_list_to_matrix
#'
#' @export
ants.image_list_to_matrix <- function(image_list, mask = NULL, sigma = NULL, epsilon = 0.5) {
  ants <- load_ants()
  ants$image_list_to_matrix(
    image_list = image_list,
    mask = mask,
    sigma = sigma,
    epsilon = epsilon
  )
}


#' @rdname ants.image_math
#'
#' @export
ants.image_math <- function(image, operation, ...) {
  ants <- load_ants()
  ants$image_math(
    image = image,
    operation = operation
  )
}


#' @rdname ants.image_mutual_information
#'
#' @export
ants.image_mutual_information <- function(image1, image2) {
  ants <- load_ants()
  ants$image_mutual_information(
    image1 = image1,
    image2 = image2
  )
}


#' @rdname ants.image_physical_space_consistency
#'
#' @export
ants.image_physical_space_consistency <- function(image1, image2, tolerance = 0.01, datatype = FALSE) {
  ants <- load_ants()
  ants$image_physical_space_consistency(
    image1 = image1,
    image2 = image2,
    tolerance = tolerance,
    datatype = datatype
  )
}


#' @rdname ants.image_read
#'
#' @export
ants.image_read <- function(filename, dimension = NULL, pixeltype = "float", reorient = FALSE) {
  ants <- load_ants()
  ants$image_read(
    filename = filename,
    dimension = dimension,
    pixeltype = pixeltype,
    reorient = reorient
  )
}


#' @rdname ants.image_similarity
#'
#' @export
ants.image_similarity <- function(fixed_image, moving_image, metric_type = "MeanSquares", fixed_mask = NULL, moving_mask = NULL, sampling_strategy = "regular", sampling_percentage = 1.0) {
  ants <- load_ants()
  ants$image_similarity(
    fixed_image = fixed_image,
    moving_image = moving_image,
    metric_type = metric_type,
    fixed_mask = fixed_mask,
    moving_mask = moving_mask,
    sampling_strategy = sampling_strategy,
    sampling_percentage = sampling_percentage
  )
}


#' @rdname ants.image_to_cluster_images
#'
#' @export
ants.image_to_cluster_images <- function(image, min_cluster_size = 50L, min_thresh = 1e-06, max_thresh = 1L) {
  ants <- load_ants()
  ants$image_to_cluster_images(
    image = image,
    min_cluster_size = min_cluster_size,
    min_thresh = min_thresh,
    max_thresh = max_thresh
  )
}


#' @rdname ants.image_type_cast
#'
#' @export
ants.image_type_cast <- function(image_list, pixeltype = NULL) {
  ants <- load_ants()
  ants$image_type_cast(
    image_list = image_list,
    pixeltype = pixeltype
  )
}


#' @rdname ants.image_write
#'
#' @export
ants.image_write <- function(image, filename, ri = FALSE) {
  ants <- load_ants()
  ants$image_write(
    image = image,
    filename = filename,
    ri = ri
  )
}


#' @rdname ants.images_from_matrix
#'
#' @export
ants.images_from_matrix <- function(data_matrix, mask) {
  ants <- load_ants()
  ants$images_from_matrix(
    data_matrix = data_matrix,
    mask = mask
  )
}


#' @rdname ants.images_to_matrix
#'
#' @export
ants.images_to_matrix <- function(image_list, mask = NULL, sigma = NULL, epsilon = 0.5) {
  ants <- load_ants()
  ants$images_to_matrix(
    image_list = image_list,
    mask = mask,
    sigma = sigma,
    epsilon = epsilon
  )
}


#' @rdname ants.iMath
#'
#' @export
ants.iMath <- function(image, operation, ...) {
  ants <- load_ants()
  ants$iMath(
    image = image,
    operation = operation
  )
}


#' @rdname ants.iMath_canny
#'
#' @export
ants.iMath_canny <- function(image, sigma, lower, upper) {
  ants <- load_ants()
  ants$iMath_canny(
    image = image,
    sigma = sigma,
    lower = lower,
    upper = upper
  )
}


#' @rdname ants.iMath_fill_holes
#'
#' @export
ants.iMath_fill_holes <- function(image, hole_type = 2L) {
  ants <- load_ants()
  ants$iMath_fill_holes(
    image = image,
    hole_type = hole_type
  )
}


#' @rdname ants.iMath_GC
#'
#' @export
ants.iMath_GC <- function(image, radius = 1L) {
  ants <- load_ants()
  ants$iMath_GC(
    image = image,
    radius = radius
  )
}


#' @rdname ants.iMath_GD
#'
#' @export
ants.iMath_GD <- function(image, radius = 1L) {
  ants <- load_ants()
  ants$iMath_GD(
    image = image,
    radius = radius
  )
}


#' @rdname ants.iMath_GE
#'
#' @export
ants.iMath_GE <- function(image, radius = 1L) {
  ants <- load_ants()
  ants$iMath_GE(
    image = image,
    radius = radius
  )
}


#' @rdname ants.iMath_get_largest_component
#'
#' @export
ants.iMath_get_largest_component <- function(image, min_size = 50L) {
  ants <- load_ants()
  ants$iMath_get_largest_component(
    image = image,
    min_size = min_size
  )
}


#' @rdname ants.iMath_GO
#'
#' @export
ants.iMath_GO <- function(image, radius = 1L) {
  ants <- load_ants()
  ants$iMath_GO(
    image = image,
    radius = radius
  )
}


#' @rdname ants.iMath_grad
#'
#' @export
ants.iMath_grad <- function(image, sigma = 0.5, normalize = FALSE) {
  ants <- load_ants()
  ants$iMath_grad(
    image = image,
    sigma = sigma,
    normalize = normalize
  )
}


#' @rdname ants.iMath_histogram_equalization
#'
#' @export
ants.iMath_histogram_equalization <- function(image, alpha, beta) {
  ants <- load_ants()
  ants$iMath_histogram_equalization(
    image = image,
    alpha = alpha,
    beta = beta
  )
}


#' @rdname ants.iMath_laplacian
#'
#' @export
ants.iMath_laplacian <- function(image, sigma = 0.5, normalize = FALSE) {
  ants <- load_ants()
  ants$iMath_laplacian(
    image = image,
    sigma = sigma,
    normalize = normalize
  )
}


#' @rdname ants.iMath_maurer_distance
#'
#' @export
ants.iMath_maurer_distance <- function(image, foreground = 1L) {
  ants <- load_ants()
  ants$iMath_maurer_distance(
    image = image,
    foreground = foreground
  )
}


#' @rdname ants.iMath_MC
#'
#' @export
ants.iMath_MC <- function(image, radius = 1L, value = 1L, shape = 1L, parametric = FALSE, lines = 3L, thickness = 1L, include_center = FALSE) {
  ants <- load_ants()
  ants$iMath_MC(
    image = image,
    radius = radius,
    value = value,
    shape = shape,
    parametric = parametric,
    lines = lines,
    thickness = thickness,
    include_center = include_center
  )
}


#' @rdname ants.iMath_MD
#'
#' @export
ants.iMath_MD <- function(image, radius = 1L, value = 1L, shape = 1L, parametric = FALSE, lines = 3L, thickness = 1L, include_center = FALSE) {
  ants <- load_ants()
  ants$iMath_MD(
    image = image,
    radius = radius,
    value = value,
    shape = shape,
    parametric = parametric,
    lines = lines,
    thickness = thickness,
    include_center = include_center
  )
}


#' @rdname ants.iMath_ME
#'
#' @export
ants.iMath_ME <- function(image, radius = 1L, value = 1L, shape = 1L, parametric = FALSE, lines = 3L, thickness = 1L, include_center = FALSE) {
  ants <- load_ants()
  ants$iMath_ME(
    image = image,
    radius = radius,
    value = value,
    shape = shape,
    parametric = parametric,
    lines = lines,
    thickness = thickness,
    include_center = include_center
  )
}


#' @rdname ants.iMath_MO
#'
#' @export
ants.iMath_MO <- function(image, radius = 1L, value = 1L, shape = 1L, parametric = FALSE, lines = 3L, thickness = 1L, include_center = FALSE) {
  ants <- load_ants()
  ants$iMath_MO(
    image = image,
    radius = radius,
    value = value,
    shape = shape,
    parametric = parametric,
    lines = lines,
    thickness = thickness,
    include_center = include_center
  )
}


#' @rdname ants.iMath_normalize
#'
#' @export
ants.iMath_normalize <- function(image) {
  ants <- load_ants()
  ants$iMath_normalize(
    image = image
  )
}


#' @rdname ants.iMath_pad
#'
#' @export
ants.iMath_pad <- function(image, padding) {
  ants <- load_ants()
  ants$iMath_pad(
    image = image,
    padding = padding
  )
}


#' @rdname ants.iMath_perona_malik
#'
#' @export
ants.iMath_perona_malik <- function(image, conductance = 0.25, n_iterations = 1L) {
  ants <- load_ants()
  ants$iMath_perona_malik(
    image = image,
    conductance = conductance,
    n_iterations = n_iterations
  )
}


#' @rdname ants.iMath_propagate_labels_through_mask
#'
#' @export
ants.iMath_propagate_labels_through_mask <- function(image, labels, stopping_value = 100L, propagation_method = 0L) {
  ants <- load_ants()
  ants$iMath_propagate_labels_through_mask(
    image = image,
    labels = labels,
    stopping_value = stopping_value,
    propagation_method = propagation_method
  )
}


#' @rdname ants.iMath_sharpen
#'
#' @export
ants.iMath_sharpen <- function(image) {
  ants <- load_ants()
  ants$iMath_sharpen(
    image = image
  )
}


#' @rdname ants.iMath_truncate_intensity
#'
#' @export
ants.iMath_truncate_intensity <- function(image, lower_q, upper_q, n_bins = 64L) {
  ants <- load_ants()
  ants$iMath_truncate_intensity(
    image = image,
    lower_q = lower_q,
    upper_q = upper_q,
    n_bins = n_bins
  )
}


#' @rdname ants.impute
#'
#' @export
ants.impute <- function(data, method = "mean", value = NULL, nan_value = nan) {
  ants <- load_ants()
  ants$impute(
    data = data,
    method = method,
    value = value,
    nan_value = nan_value
  )
}


#' @rdname ants.initialize_eigenanatomy
#'
#' @export
ants.initialize_eigenanatomy <- function(initmat, mask = NULL, initlabels = NULL, nreps = 1L, smoothing = 0L) {
  ants <- load_ants()
  ants$initialize_eigenanatomy(
    initmat = initmat,
    mask = mask,
    initlabels = initlabels,
    nreps = nreps,
    smoothing = smoothing
  )
}


#' @rdname ants.integrate_velocity_field
#'
#' @export
ants.integrate_velocity_field <- function(velocity_field, lower_integration_bound = 0.0, upper_integration_bound = 1.0, number_of_integration_steps = 10L) {
  ants <- load_ants()
  ants$integrate_velocity_field(
    velocity_field = velocity_field,
    lower_integration_bound = lower_integration_bound,
    upper_integration_bound = upper_integration_bound,
    number_of_integration_steps = number_of_integration_steps
  )
}


#' @rdname ants.invert_ants_transform
#'
#' @export
ants.invert_ants_transform <- function(transform) {
  ants <- load_ants()
  ants$invert_ants_transform(
    transform = transform
  )
}


#' @rdname ants.invert_displacement_field
#'
#' @export
ants.invert_displacement_field <- function(displacement_field, inverse_field_initial_estimate, maximum_number_of_iterations = 20L, mean_error_tolerance_threshold = 0.001, max_error_tolerance_threshold = 0.1, enforce_boundary_condition = TRUE) {
  ants <- load_ants()
  ants$invert_displacement_field(
    displacement_field = displacement_field,
    inverse_field_initial_estimate = inverse_field_initial_estimate,
    maximum_number_of_iterations = maximum_number_of_iterations,
    mean_error_tolerance_threshold = mean_error_tolerance_threshold,
    max_error_tolerance_threshold = max_error_tolerance_threshold,
    enforce_boundary_condition = enforce_boundary_condition
  )
}


#' @rdname ants.joint_label_fusion
#'
#' @export
ants.joint_label_fusion <- function(target_image, target_image_mask, atlas_list, beta = 4L, rad = 2L, label_list = NULL, rho = 0.01, usecor = FALSE, r_search = 3L, nonnegative = FALSE, no_zeroes = FALSE, max_lab_plus_one = FALSE, output_prefix = NULL, verbose = FALSE) {
  ants <- load_ants()
  ants$joint_label_fusion(
    target_image = target_image,
    target_image_mask = target_image_mask,
    atlas_list = atlas_list,
    beta = beta,
    rad = rad,
    label_list = label_list,
    rho = rho,
    usecor = usecor,
    r_search = r_search,
    nonnegative = nonnegative,
    no_zeroes = no_zeroes,
    max_lab_plus_one = max_lab_plus_one,
    output_prefix = output_prefix,
    verbose = verbose
  )
}


#' @rdname ants.kelly_kapowski
#'
#' @export
ants.kelly_kapowski <- function(s, g, w, its = 45L, r = 0.025, m = 1.5) {
  ants <- load_ants()
  ants$kelly_kapowski(
    s = s,
    g = g,
    w = w,
    its = its,
    r = r,
    m = m
  )
}


#' @rdname ants.kmeans_segmentation
#'
#' @export
ants.kmeans_segmentation <- function(image, k, kmask = NULL, mrf = 0.1) {
  ants <- load_ants()
  ants$kmeans_segmentation(
    image = image,
    k = k,
    kmask = kmask,
    mrf = mrf
  )
}


#' @rdname ants.label_clusters
#'
#' @export
ants.label_clusters <- function(image, min_cluster_size = 50L, min_thresh = 1e-06, max_thresh = 1L, fully_connected = FALSE) {
  ants <- load_ants()
  ants$label_clusters(
    image = image,
    min_cluster_size = min_cluster_size,
    min_thresh = min_thresh,
    max_thresh = max_thresh,
    fully_connected = fully_connected
  )
}


#' @rdname ants.label_geometry_measures
#'
#' @export
ants.label_geometry_measures <- function(label_image, intensity_image = NULL) {
  ants <- load_ants()
  ants$label_geometry_measures(
    label_image = label_image,
    intensity_image = intensity_image
  )
}


#' @rdname ants.label_image_centroids
#'
#' @export
ants.label_image_centroids <- function(image, physical = FALSE, convex = TRUE, verbose = FALSE) {
  ants <- load_ants()
  ants$label_image_centroids(
    image = image,
    physical = physical,
    convex = convex,
    verbose = verbose
  )
}


#' @rdname ants.label_overlap_measures
#'
#' @export
ants.label_overlap_measures <- function(source_image, target_image) {
  ants <- load_ants()
  ants$label_overlap_measures(
    source_image = source_image,
    target_image = target_image
  )
}


#' @rdname ants.label_stats
#'
#' @export
ants.label_stats <- function(image, label_image) {
  ants <- load_ants()
  ants$label_stats(
    image = image,
    label_image = label_image
  )
}


#' @rdname ants.LabelImage
#'
#' @export
ants.LabelImage <- function(label_image, label_info = NULL, template = NULL) {
  ants <- load_ants()
  ants$LabelImage(
    label_image = label_image,
    label_info = label_info,
    template = template
  )
}



#' @rdname ants.list_to_ndimage
#'
#' @export
ants.list_to_ndimage <- function(image, image_list) {
  ants <- load_ants()
  ants$list_to_ndimage(
    image = image,
    image_list = image_list
  )
}


#' @rdname ants.local_joint_label_fusion
#'
#' @export
ants.local_joint_label_fusion <- function(target_image, which_labels, target_mask, initial_label, atlas_list, label_list, submask_dilation = 10L, type_of_transform = "SyN", aff_metric = "meansquares", syn_metric = "mattes", syn_sampling = 32L, reg_iterations = list(40L, 20L, 0L), aff_iterations = list(500L, 50L, 0L), grad_step = 0.2, flow_sigma = 3L, total_sigma = 0L, beta = 4L, rad = 2L, rho = 0.1, usecor = FALSE, r_search = 3L, nonnegative = FALSE, no_zeroes = FALSE, max_lab_plus_one = FALSE, local_mask_transform = "Similarity", output_prefix = NULL, verbose = FALSE) {
  ants <- load_ants()
  ants$local_joint_label_fusion(
    target_image = target_image,
    which_labels = which_labels,
    target_mask = target_mask,
    initial_label = initial_label,
    atlas_list = atlas_list,
    label_list = label_list,
    submask_dilation = submask_dilation,
    type_of_transform = type_of_transform,
    aff_metric = aff_metric,
    syn_metric = syn_metric,
    syn_sampling = syn_sampling,
    reg_iterations = reg_iterations,
    aff_iterations = aff_iterations,
    grad_step = grad_step,
    flow_sigma = flow_sigma,
    total_sigma = total_sigma,
    beta = beta,
    rad = rad,
    rho = rho,
    usecor = usecor,
    r_search = r_search,
    nonnegative = nonnegative,
    no_zeroes = no_zeroes,
    max_lab_plus_one = max_lab_plus_one,
    local_mask_transform = local_mask_transform,
    output_prefix = output_prefix,
    verbose = verbose
  )
}


#' @rdname ants.make_image
#'
#' @export
ants.make_image <- function(imagesize, voxval = 0L, spacing = NULL, origin = NULL, direction = NULL, has_components = FALSE, pixeltype = "float") {
  ants <- load_ants()
  ants$make_image(
    imagesize = imagesize,
    voxval = voxval,
    spacing = spacing,
    origin = origin,
    direction = direction,
    has_components = has_components,
    pixeltype = pixeltype
  )
}


#' @rdname ants.make_points_image
#'
#' @export
ants.make_points_image <- function(pts, mask, radius = 5L) {
  ants <- load_ants()
  ants$make_points_image(
    pts = pts,
    mask = mask,
    radius = radius
  )
}


#' @rdname ants.mask_image
#'
#' @export
ants.mask_image <- function(image, mask, level = 1L, binarize = FALSE) {
  ants <- load_ants()
  ants$mask_image(
    image = image,
    mask = mask,
    level = level,
    binarize = binarize
  )
}


#' @rdname ants.matrix_from_images
#'
#' @export
ants.matrix_from_images <- function(image_list, mask = NULL, sigma = NULL, epsilon = 0.5) {
  ants <- load_ants()
  ants$matrix_from_images(
    image_list = image_list,
    mask = mask,
    sigma = sigma,
    epsilon = epsilon
  )
}


#' @rdname ants.matrix_to_images
#'
#' @export
ants.matrix_to_images <- function(data_matrix, mask) {
  ants <- load_ants()
  ants$matrix_to_images(
    data_matrix = data_matrix,
    mask = mask
  )
}


#' @rdname ants.matrix_to_timeseries
#'
#' @export
ants.matrix_to_timeseries <- function(image, matrix, mask = NULL) {
  ants <- load_ants()
  ants$matrix_to_timeseries(
    image = image,
    matrix = matrix,
    mask = mask
  )
}


#' @rdname ants.merge_channels
#'
#' @export
ants.merge_channels <- function(image_list) {
  ants <- load_ants()
  ants$merge_channels(
    image_list = image_list
  )
}


#' @rdname ants.mni2tal
#'
#' @export
ants.mni2tal <- function(xin) {
  ants <- load_ants()
  ants$mni2tal(
    xin = xin
  )
}


#' @rdname ants.morphology
#'
#' @export
ants.morphology <- function(image, operation, radius, mtype = "binary", value = 1L, shape = "ball", radius_is_parametric = FALSE, thickness = 1L, lines = 3L, include_center = FALSE) {
  ants <- load_ants()
  ants$morphology(
    image = image,
    operation = operation,
    radius = radius,
    mtype = mtype,
    value = value,
    shape = shape,
    radius_is_parametric = radius_is_parametric,
    thickness = thickness,
    lines = lines,
    include_center = include_center
  )
}


#' @rdname ants.motion_correction
#'
#' @export
ants.motion_correction <- function(image, fixed = NULL, type_of_transform = "BOLDRigid", mask = NULL, fdOffset = 50L, outprefix = "", verbose = FALSE) {
  ants <- load_ants()
  ants$motion_correction(
    image = image,
    fixed = fixed,
    type_of_transform = type_of_transform,
    mask = mask,
    fdOffset = fdOffset,
    outprefix = outprefix,
    verbose = verbose
  )
}


#' @rdname ants.movie
#'
#' @export
ants.movie <- function(image, filename = NULL, writer = NULL, fps = 30L) {
  ants <- load_ants()
  ants$movie(
    image = image,
    filename = filename,
    writer = writer,
    fps = fps
  )
}


#' @rdname ants.multi_label_morphology
#'
#' @export
ants.multi_label_morphology <- function(image, operation, radius, dilation_mask = NULL, label_list = NULL, force = FALSE) {
  ants <- load_ants()
  ants$multi_label_morphology(
    image = image,
    operation = operation,
    radius = radius,
    dilation_mask = dilation_mask,
    label_list = label_list,
    force = force
  )
}


#' @rdname ants.multiply_images
#'
#' @export
ants.multiply_images <- function(image1, image2) {
  ants <- load_ants()
  ants$multiply_images(
    image1 = image1,
    image2 = image2
  )
}


#' @rdname ants.n3_bias_field_correction
#'
#' @export
ants.n3_bias_field_correction <- function(image, downsample_factor = 3L) {
  ants <- load_ants()
  ants$n3_bias_field_correction(
    image = image,
    downsample_factor = downsample_factor
  )
}


#' @rdname ants.ndimage_to_list
#'
#' @export
ants.ndimage_to_list <- function(image) {
  ants <- load_ants()
  ants$ndimage_to_list(
    image = image
  )
}


#' @rdname ants.new_ants_metric
#'
#' @export
ants.new_ants_metric <- function(dimension = 3L, precision = "float", metric_type = "MeanSquares") {
  ants <- load_ants()
  ants$new_ants_metric(
    dimension = dimension,
    precision = precision,
    metric_type = metric_type
  )
}


#' @rdname ants.new_ants_transform
#'
#' @export
ants.new_ants_transform <- function(precision = "float", dimension = 3L, transform_type = "AffineTransform", parameters = NULL) {
  ants <- load_ants()
  ants$new_ants_transform(
    precision = precision,
    dimension = dimension,
    transform_type = transform_type,
    parameters = parameters
  )
}


#' @rdname ants.nifti_to_ants
#'
#' @export
ants.nifti_to_ants <- function(nib_image) {
  ants <- load_ants()
  ants$nifti_to_ants(
    nib_image = nib_image
  )
}


#' @rdname ants.otsu_segmentation
#'
#' @export
ants.otsu_segmentation <- function(image, k, mask = NULL) {
  ants <- load_ants()
  ants$otsu_segmentation(
    image = image,
    k = k,
    mask = mask
  )
}


#' @rdname ants.pad_image
#'
#' @export
ants.pad_image <- function(image, shape = NULL, pad_width = NULL, value = 0.0, return_padvals = FALSE) {
  ants <- load_ants()
  ants$pad_image(
    image = image,
    shape = shape,
    pad_width = pad_width,
    value = value,
    return_padvals = return_padvals
  )
}





#' @rdname ants.prior_based_segmentation
#'
#' @export
ants.prior_based_segmentation <- function(image, priors, mask, priorweight = 0.25, mrf = 0.1, iterations = 25L) {
  ants <- load_ants()
  ants$prior_based_segmentation(
    image = image,
    priors = priors,
    mask = mask,
    priorweight = priorweight,
    mrf = mrf,
    iterations = iterations
  )
}


#' @rdname ants.quantile
#'
#' @export
ants.quantile <- function(image, q, nonzero = TRUE) {
  ants <- load_ants()
  ants$quantile(
    image = image,
    q = q,
    nonzero = nonzero
  )
}


#' @rdname ants.rank_intensity
#'
#' @export
ants.rank_intensity <- function(x, mask = NULL, get_mask = TRUE, method = "max") {
  ants <- load_ants()
  ants$rank_intensity(
    x = x,
    mask = mask,
    get_mask = get_mask,
    method = method
  )
}


#' @rdname ants.read_transform
#'
#' @export
ants.read_transform <- function(filename, precision = "float") {
  ants <- load_ants()
  ants$read_transform(
    filename = filename,
    precision = precision
  )
}


#' @rdname ants.reflect_image
#'
#' @export
ants.reflect_image <- function(image, axis = NULL, tx = NULL, metric = "mattes") {
  ants <- load_ants()
  ants$reflect_image(
    image = image,
    axis = axis,
    tx = tx,
    metric = metric
  )
}


#' @rdname ants.registration
#'
#' @export
ants.registration <- function(fixed, moving, type_of_transform = "SyN", initial_transform = NULL, outprefix = "", mask = NULL, grad_step = 0.2, flow_sigma = 3L, total_sigma = 0L, aff_metric = "mattes", aff_sampling = 32L, aff_random_sampling_rate = 0.2, syn_metric = "mattes", syn_sampling = 32L, reg_iterations = list(40L, 20L, 0L), aff_iterations = list(2100L, 1200L, 1200L, 10L), aff_shrink_factors = list(6L, 4L, 2L, 1L), aff_smoothing_sigmas = list(3L, 2L, 1L, 0L), write_composite_transform = FALSE, random_seed = NULL, verbose = FALSE, multivariate_extras = NULL, restrict_transformation = NULL, smoothing_in_mm = FALSE) {
  ants <- load_ants()
  ants$registration(
    fixed = fixed,
    moving = moving,
    type_of_transform = type_of_transform,
    initial_transform = initial_transform,
    outprefix = outprefix,
    mask = mask,
    grad_step = grad_step,
    flow_sigma = flow_sigma,
    total_sigma = total_sigma,
    aff_metric = aff_metric,
    aff_sampling = aff_sampling,
    aff_random_sampling_rate = aff_random_sampling_rate,
    syn_metric = syn_metric,
    syn_sampling = syn_sampling,
    reg_iterations = reg_iterations,
    aff_iterations = aff_iterations,
    aff_shrink_factors = aff_shrink_factors,
    aff_smoothing_sigmas = aff_smoothing_sigmas,
    write_composite_transform = write_composite_transform,
    random_seed = random_seed,
    verbose = verbose,
    multivariate_extras = multivariate_extras,
    restrict_transformation = restrict_transformation,
    smoothing_in_mm = smoothing_in_mm
  )
}


#' @rdname ants.regress_components
#'
#' @export
ants.regress_components <- function(data, components, remove_mean = TRUE) {
  ants <- load_ants()
  ants$regress_components(
    data = data,
    components = components,
    remove_mean = remove_mean
  )
}


#' @rdname ants.regress_poly
#'
#' @export
ants.regress_poly <- function(degree, data, remove_mean = TRUE, axis = -1L) {
  ants <- load_ants()
  ants$regress_poly(
    degree = degree,
    data = data,
    remove_mean = remove_mean,
    axis = axis
  )
}


#' @rdname ants.render_surface_function
#'
#' @export
ants.render_surface_function <- function(surfimg, funcimg = NULL, alphasurf = 0.2, alphafunc = 1.0, isosurf = 0.5, isofunc = 0.5, smoothsurf = NULL, smoothfunc = NULL, cmapsurf = "grey", cmapfunc = "red", filename = NULL, notebook = FALSE, auto_open = FALSE) {
  ants <- load_ants()
  ants$render_surface_function(
    surfimg = surfimg,
    funcimg = funcimg,
    alphasurf = alphasurf,
    alphafunc = alphafunc,
    isosurf = isosurf,
    isofunc = isofunc,
    smoothsurf = smoothsurf,
    smoothfunc = smoothfunc,
    cmapsurf = cmapsurf,
    cmapfunc = cmapfunc,
    filename = filename,
    notebook = notebook,
    auto_open = auto_open
  )
}


#' @rdname ants.reorient_image2
#'
#' @export
ants.reorient_image2 <- function(image, orientation = "RAS") {
  ants <- load_ants()
  ants$reorient_image2(
    image = image,
    orientation = orientation
  )
}


#' @rdname ants.resample_image
#'
#' @export
ants.resample_image <- function(image, resample_params, use_voxels = FALSE, interp_type = 1L) {
  ants <- load_ants()
  ants$resample_image(
    image = image,
    resample_params = resample_params,
    use_voxels = use_voxels,
    interp_type = interp_type
  )
}


#' @rdname ants.resample_image_to_target
#'
#' @export
ants.resample_image_to_target <- function(image, target, interp_type = "linear", imagetype = 0L, verbose = FALSE) {
  ants <- load_ants()
  ants$resample_image_to_target(
    image = image,
    target = target,
    interp_type = interp_type,
    imagetype = imagetype,
    verbose = verbose
  )
}


#' @rdname ants.rgb_to_vector
#'
#' @export
ants.rgb_to_vector <- function(image) {
  ants <- load_ants()
  ants$rgb_to_vector(
    image = image
  )
}


#' @rdname ants.scalar_to_rgb
#'
#' @export
ants.scalar_to_rgb <- function(image, mask = NULL, filename = NULL, cmap = "red", custom_colormap_file = NULL, min_input = NULL, max_input = NULL, min_rgb_output = NULL, max_rgb_output = NULL, vtk_lookup_table = NULL) {
  ants <- load_ants()
  ants$scalar_to_rgb(
    image = image,
    mask = mask,
    filename = filename,
    cmap = cmap,
    custom_colormap_file = custom_colormap_file,
    min_input = min_input,
    max_input = max_input,
    min_rgb_output = min_rgb_output,
    max_rgb_output = max_rgb_output,
    vtk_lookup_table = vtk_lookup_table
  )
}


#' @rdname ants.set_ants_transform_fixed_parameters
#'
#' @export
ants.set_ants_transform_fixed_parameters <- function(transform, parameters) {
  ants <- load_ants()
  ants$set_ants_transform_fixed_parameters(
    transform = transform,
    parameters = parameters
  )
}


#' @rdname ants.set_ants_transform_parameters
#'
#' @export
ants.set_ants_transform_parameters <- function(transform, parameters) {
  ants <- load_ants()
  ants$set_ants_transform_parameters(
    transform = transform,
    parameters = parameters
  )
}


#' @rdname ants.set_direction
#'
#' @export
ants.set_direction <- function(image, direction) {
  ants <- load_ants()
  ants$set_direction(
    image = image,
    direction = direction
  )
}


#' @rdname ants.set_origin
#'
#' @export
ants.set_origin <- function(image, origin) {
  ants <- load_ants()
  ants$set_origin(
    image = image,
    origin = origin
  )
}


#' @rdname ants.set_spacing
#'
#' @export
ants.set_spacing <- function(image, spacing) {
  ants <- load_ants()
  ants$set_spacing(
    image = image,
    spacing = spacing
  )
}


#' @rdname ants.short_ptype
#'
#' @export
ants.short_ptype <- function(pixeltype) {
  ants <- load_ants()
  ants$short_ptype(
    pixeltype = pixeltype
  )
}


#' @rdname ants.simulate_displacement_field
#'
#' @export
ants.simulate_displacement_field <- function(domain_image, field_type = "bspline", number_of_random_points = 1000L, sd_noise = 10.0, enforce_stationary_boundary = TRUE, number_of_fitting_levels = 4L, mesh_size = 1L, sd_smoothing = 4.0) {
  ants <- load_ants()
  ants$simulate_displacement_field(
    domain_image = domain_image,
    field_type = field_type,
    number_of_random_points = number_of_random_points,
    sd_noise = sd_noise,
    enforce_stationary_boundary = enforce_stationary_boundary,
    number_of_fitting_levels = number_of_fitting_levels,
    mesh_size = mesh_size,
    sd_smoothing = sd_smoothing
  )
}


#' @rdname ants.slice_image
#'
#' @export
ants.slice_image <- function(image, axis = NULL, idx = NULL, collapse_strategy = 0L) {
  ants <- load_ants()
  ants$slice_image(
    image = image,
    axis = axis,
    idx = idx,
    collapse_strategy = collapse_strategy
  )
}


#' @rdname ants.smooth_image
#'
#' @export
ants.smooth_image <- function(image, sigma, sigma_in_physical_coordinates = TRUE, FWHM = FALSE, max_kernel_width = 32L) {
  ants <- load_ants()
  ants$smooth_image(
    image = image,
    sigma = sigma,
    sigma_in_physical_coordinates = sigma_in_physical_coordinates,
    FWHM = FWHM,
    max_kernel_width = max_kernel_width
  )
}



#' @rdname ants.split_channels
#'
#' @export
ants.split_channels <- function(image) {
  ants <- load_ants()
  ants$split_channels(
    image = image
  )
}


#' @rdname ants.supported_metrics
#'
#' @export
ants.supported_metrics <- function() {
  ants <- load_ants()
  ants$supported_metrics(
)
}


#' @rdname ants.surf
#'
#' @export
ants.surf <- function(x, y = NULL, z = NULL, quantlimits = list(0.1, 0.9), colormap = "jet", grayscale = 0.7, bg_grayscale = 0.9, alpha = NULL, inflation_factor = 0L, tol = 0.03, smoothing_sigma = 0.0, rotation_params = list(90L, 0L, 270L), overlay_limits = NULL, filename = NULL, verbose = FALSE) {
  ants <- load_ants()
  ants$surf(
    x = x,
    y = y,
    z = z,
    quantlimits = quantlimits,
    colormap = colormap,
    grayscale = grayscale,
    bg_grayscale = bg_grayscale,
    alpha = alpha,
    inflation_factor = inflation_factor,
    tol = tol,
    smoothing_sigma = smoothing_sigma,
    rotation_params = rotation_params,
    overlay_limits = overlay_limits,
    filename = filename,
    verbose = verbose
  )
}


#' @rdname ants.surf_fold
#'
#' @export
ants.surf_fold <- function(image, outfile, dilation = 0L, inflation = 10L, alpha = 1.0, overlay = NULL, overlay_mask = NULL, overlay_cmap = "jet", overlay_scale = FALSE, overlay_alpha = 1.0, rotation = NULL, cut_idx = NULL, cut_side = "left", grayscale = 0.7, bg_grayscale = 0.9, verbose = FALSE, cleanup = TRUE) {
  ants <- load_ants()
  ants$surf_fold(
    image = image,
    outfile = outfile,
    dilation = dilation,
    inflation = inflation,
    alpha = alpha,
    overlay = overlay,
    overlay_mask = overlay_mask,
    overlay_cmap = overlay_cmap,
    overlay_scale = overlay_scale,
    overlay_alpha = overlay_alpha,
    rotation = rotation,
    cut_idx = cut_idx,
    cut_side = cut_side,
    grayscale = grayscale,
    bg_grayscale = bg_grayscale,
    verbose = verbose,
    cleanup = cleanup
  )
}


#' @rdname ants.surf_smooth
#'
#' @export
ants.surf_smooth <- function(image, outfile, dilation = 1.0, smooth = 1.0, threshold = 0.5, inflation = 200L, alpha = 1.0, cut_idx = NULL, cut_side = "left", overlay = NULL, overlay_mask = NULL, overlay_cmap = "jet", overlay_scale = FALSE, overlay_alpha = 1.0, rotation = NULL, grayscale = 0.7, bg_grayscale = 0.9, verbose = FALSE, cleanup = TRUE) {
  ants <- load_ants()
  ants$surf_smooth(
    image = image,
    outfile = outfile,
    dilation = dilation,
    smooth = smooth,
    threshold = threshold,
    inflation = inflation,
    alpha = alpha,
    cut_idx = cut_idx,
    cut_side = cut_side,
    overlay = overlay,
    overlay_mask = overlay_mask,
    overlay_cmap = overlay_cmap,
    overlay_scale = overlay_scale,
    overlay_alpha = overlay_alpha,
    rotation = rotation,
    grayscale = grayscale,
    bg_grayscale = bg_grayscale,
    verbose = verbose,
    cleanup = cleanup
  )
}


#' @rdname ants.symmetrize_image
#'
#' @export
ants.symmetrize_image <- function(image) {
  ants <- load_ants()
  ants$symmetrize_image(
    image = image
  )
}


#' @rdname ants.threshold_image
#'
#' @export
ants.threshold_image <- function(image, low_thresh = NULL, high_thresh = NULL, inval = 1L, outval = 0L, binary = TRUE) {
  ants <- load_ants()
  ants$threshold_image(
    image = image,
    low_thresh = low_thresh,
    high_thresh = high_thresh,
    inval = inval,
    outval = outval,
    binary = binary
  )
}


#' @rdname ants.timeseries_to_matrix
#'
#' @export
ants.timeseries_to_matrix <- function(image, mask = NULL) {
  ants <- load_ants()
  ants$timeseries_to_matrix(
    image = image,
    mask = mask
  )
}


#' @rdname ants.to_nibabel
#'
#' @export
ants.to_nibabel <- function(image) {
  ants <- load_ants()
  ants$to_nibabel(
    image = image
  )
}


#' @rdname ants.transform_from_displacement_field
#'
#' @export
ants.transform_from_displacement_field <- function(field) {
  ants <- load_ants()
  ants$transform_from_displacement_field(
    field = field
  )
}


#' @rdname ants.transform_index_to_physical_point
#'
#' @export
ants.transform_index_to_physical_point <- function(image, index) {
  ants <- load_ants()
  ants$transform_index_to_physical_point(
    image = image,
    index = index
  )
}


#' @rdname ants.transform_physical_point_to_index
#'
#' @export
ants.transform_physical_point_to_index <- function(image, point) {
  ants <- load_ants()
  ants$transform_physical_point_to_index(
    image = image,
    point = point
  )
}


#' @rdname ants.transform_to_displacement_field
#'
#' @export
ants.transform_to_displacement_field <- function(xfrm, ref) {
  ants <- load_ants()
  ants$transform_to_displacement_field(
    xfrm = xfrm,
    ref = ref
  )
}


#' @rdname ants.vector_to_rgb
#'
#' @export
ants.vector_to_rgb <- function(image) {
  ants <- load_ants()
  ants$vector_to_rgb(
    image = image
  )
}


#' @rdname ants.vol
#'
#' @export
ants.vol <- function(volume, overlays = NULL, quantlimits = list(0.1, 0.9), colormap = "jet", rotation_params = list(90L, 0L, 270L), overlay_limits = NULL, magnification_factor = 1.0, intensity_truncation = list(0.0, 1.0), filename = NULL, verbose = FALSE) {
  ants <- load_ants()
  ants$vol(
    volume = volume,
    overlays = overlays,
    quantlimits = quantlimits,
    colormap = colormap,
    rotation_params = rotation_params,
    overlay_limits = overlay_limits,
    magnification_factor = magnification_factor,
    intensity_truncation = intensity_truncation,
    filename = filename,
    verbose = verbose
  )
}


#' @rdname ants.vol_fold
#'
#' @export
ants.vol_fold <- function(image, outfile, magnification = 1.0, dilation = 0L, inflation = 10L, alpha = 1.0, overlay = NULL, overlay_mask = NULL, overlay_cmap = "jet", overlay_scale = FALSE, overlay_alpha = 1.0, rotation = NULL, cut_idx = NULL, cut_side = "left", grayscale = 0.7, bg_grayscale = 0.9, verbose = FALSE, cleanup = TRUE) {
  ants <- load_ants()
  ants$vol_fold(
    image = image,
    outfile = outfile,
    magnification = magnification,
    dilation = dilation,
    inflation = inflation,
    alpha = alpha,
    overlay = overlay,
    overlay_mask = overlay_mask,
    overlay_cmap = overlay_cmap,
    overlay_scale = overlay_scale,
    overlay_alpha = overlay_alpha,
    rotation = rotation,
    cut_idx = cut_idx,
    cut_side = cut_side,
    grayscale = grayscale,
    bg_grayscale = bg_grayscale,
    verbose = verbose,
    cleanup = cleanup
  )
}


#' @rdname ants.weingarten_image_curvature
#'
#' @export
ants.weingarten_image_curvature <- function(image, sigma = 1.0, opt = "mean") {
  ants <- load_ants()
  ants$weingarten_image_curvature(
    image = image,
    sigma = sigma,
    opt = opt
  )
}


#' @rdname ants.write_transform
#'
#' @export
ants.write_transform <- function(transform, filename) {
  ants <- load_ants()
  ants$write_transform(
    transform = transform,
    filename = filename
  )
}

