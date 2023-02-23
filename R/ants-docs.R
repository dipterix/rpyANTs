
#' @name ants.abp_n4
#'
#' @title Truncate outlier intensities and bias correct with the 'N4' algorithm.
#' @param image 'ANTsImage' image to correct and truncate
#' @param intensity_truncation 3-list quantiles for intensity truncation
#' @param mask 'ANTsImage' (optional) mask for bias correction
#' @param usen3 logical, if true, use 'N3' bias correction instead of 'N4' Returns
#' @examples
#' if( interactive() && ants_available() ) {
#'
#' image <- ants.image_read(ants.get_ants_data('r16'))
#' image2 = ants.abp_n4(image)
#'
#' ants.plot(image)
#' ants.plot(image2)
#' }
NULL



#' @name ants.add_noise_to_image
#' @title Add noise to an image
#' @description
#' Add noise to an image using additive Guassian, salt-and-pepper, shot, or speckle noise.
#'
#' @param image 'ANTsImage' scalar image
#' @param noise_model string \code{'additivegaussian'}, \code{'saltandpepper'},
#' \code{'shot'}, or \code{'speckle'}
#' @param noise_parameters \code{\link{tuple}} or array or float, depending on
#' \code{noise_model}. For \code{'additivegaussian'}, this is (mean,
#' standard deviation), for \code{'saltandpepper'} (probability, salt-value,
#' pepper-value), for \code{'shot'} (scale), for \code{'speckle'} (standard
#' deviation)
#'
#' @examples
#'
#' if( interactive() && ants_available() ) {
#'
#' image = ants.image_read(ants.get_ants_data('r16'))
#' noise_image1 = ants.add_noise_to_image(
#'   image, 'additivegaussian', list(0.0, 1.0))
#' noise_image2 = ants.add_noise_to_image(
#'   image, 'saltandpepper', list(0.1, 0.0, 100.0))
#' noise_image3 = ants.add_noise_to_image(image, 'shot', 1.0)
#' noise_image4 = ants.add_noise_to_image(image, 'speckle', 1.0)
#'
#' }
#'
NULL




#' @name ants.affine_initializer
#' @title A multi-start optimizer for affine registration
#' @description
#' Searches over the sphere to find a good initialization for further
#' registration refinement, if needed. This is a arapper for the ANTs
#' function antsAffineInitializer. ANTsR function: `affineInitializer` Arguments
#'
#' @param fixed_image : ANTsImage:
#' the fixed reference image moving_image : ANTsImage the moving image to be mapped to the fixed space search_factor : scalar degree of increments on the sphere to search radian_fraction : scalar between zero and one, defines the arc to search over use_principal_axis : boolean boolean to initialize by principal axis local_search_iterations : scalar gradient descent iterations mask : ANTsImage (optional) optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param moving_image : ANTsImage :
#' the moving image to be mapped to the fixed space search_factor : scalar degree of increments on the sphere to search radian_fraction : scalar between zero and one, defines the arc to search over use_principal_axis : boolean boolean to initialize by principal axis local_search_iterations : scalar gradient descent iterations mask : ANTsImage (optional) optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param search_factor : scalar:
#' degree of increments on the sphere to search radian_fraction : scalar between zero and one, defines the arc to search over use_principal_axis : boolean boolean to initialize by principal axis local_search_iterations : scalar gradient descent iterations mask : ANTsImage (optional) optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param radian_fraction : scalar:
#' between zero and one, defines the arc to search over use_principal_axis : boolean boolean to initialize by principal axis local_search_iterations : scalar gradient descent iterations mask : ANTsImage (optional) optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param use_principal_axis : boolean:
#' boolean to initialize by principal axis local_search_iterations : scalar gradient descent iterations mask : ANTsImage (optional) optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param local_search_iterations : scalar:
#' gradient descent iterations mask : ANTsImage (optional) optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param mask : ANTsImage (optional):
#' optional mask to restrict registration txfn : string (optional) filename for the transformation
#'
#' @param txfn : string (optional):
#' filename for the transformation
#'
#' @examples
#'
#' if( interactive() && ants_available() ) {
#'
#' fi <- ants.image_read(ants.get_ants_data('r16'))
#' mi <- ants.image_read(ants.get_ants_data('r27'))
#' txfile <- ants.affine_initializer( fi, mi )
#' tx <- ants.read_transform(txfile)
#' tx
#'
#' }
#'
NULL


#' @name ants.allclose
#' @title Check if two images have the same array values
#' @param image1 image1
#' @param image2 image2
NULL



#' @name ants.anti_alias
#'
#' @title Apply Anti-Alias filter to a binary image
#' @param image : ANTsImage:
#' binary image to which anti-aliasing will be applied
#'
#' @examples
#'
#' if( interactive() && ants_available() ) {
#'
#'
#' img <- ants.image_read(ants.get_data('r16'))
#' mask <- ants.get_mask(img)
#' mask_aa <- ants.anti_alias(img)
#' mask_aa
#'
#'
#' }
#'
NULL


#' @name ants.ANTsImage
#'
#' @title ANTs Class ANTsImage
#' @param pixeltype pixeltype
#' @param dimension dimension
#' @param components components
#' @param pointer pointer
#' @param is_rgb is_rgb
#' @param label_image label_image
NULL



#' @name ants.ANTsTransform
#'
#' @title ANTs Class ANTsTransform
#' @param precision precision
#' @param dimension dimension
#' @param transform_type transform_type
#' @param pointer pointer
NULL


#' @name ants.apply_ants_transform
#'
#' @title Apply ANTsTransform to data
#' @param transform : ANTsTransform:
#' transform to apply to image
#'
#' @param data : ndarray/list/tuple:
#' data to which transform will be applied
#'
#' @param data_type : string:
#' type of data Options : 'point' 'vector' 'image'
#'
#' @param reference : ANTsImage:
#' target space for transforming image
#'
#' @param interpolation : string:
#' type of interpolation to use
#'
#' @param point : list/tuple:
#' point to which the transform will be applied
#'
#' @param vector : list/tuple:
#' vector to which the transform will be applied
#'
#' @param ... additional options passed to `apply_ants_transform_to_image`
#'
#' @examples
#'
#' if( interactive() && ants_available() ) {
#'
#' img <- ants.image_read(ants.get_ants_data("r16"))$clone('float')
#' tx <- ants.new_ants_transform(dimension=2L)
#' tx$set_parameters(tuple(0.9,0,0,1.1,10,11))
#' img2 <- tx$apply_to_image(img, img)
#'
#' ants.plot(img)
#' ants.plot(img2)
#'
#'
#' tx <- ants.new_ants_transform()
#' params <- tx$parameters
#' tx$set_parameters(params*2)
#' pt2 <- tx$apply_to_point(tuple(1,2,3))
#'
#' # should be (2,4,6)
#' np <- import("numpy")
#' np$asarray(pt2)
#'
#'
#' }
#'
NULL

