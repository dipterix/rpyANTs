#' @name ants.plot
#' @title Plot an 'ANTsImage'
#' @description
#' Use mask_image and/or threshold_image to preprocess images to be be
#' overlayed and display the overlays in a given range. See the wiki examples.
#' By default, images will be reoriented to 'LAI' orientation before plotting.
#' So, if axis == 0, the images will be ordered from the
#' left side of the brain to the right side of the brain. If axis == 1,
#' the images will be ordered from the anterior (front) of the brain to
#' the posterior (back) of the brain. And if axis == 2, the images will
#' be ordered from the inferior (bottom) of the brain to the superior (top)
#' of the brain.
#'
#' @param image 'ANTsImage', image to plot
#' @param overlay 'ANTsImage', image to overlay
#' @param blend blend
#' @param alpha alpha
#' @param cmap string, colormap to use for base image. See matplotlib.
#' @param overlay_cmap string, colormap to use for overlay images, if applicable. See matplotlib.
#' @param overlay_alpha float, level of transparency for any overlays. Smaller value means the overlay is more transparent. See matplotlib.
#' @param vminol vminol
#' @param vmaxol vmaxol
#' @param cbar cbar
#' @param cbar_length cbar_length
#' @param cbar_dx cbar_dx
#' @param cbar_vertical cbar_vertical
#' @param axis integer, which axis to plot along if image is 3D
#' @param nslices integer, number of slices to plot if image is 3D
#' @param slices slices : list or \code{\link{tuple}} of integers, specific
#' slice indices to plot if image is 3D. If given, this will override
#' \code{nslices}. This can be absolute array indices (e.g.
#' \code{list(80,100,120)}), or this can be relative array indices
#' (e.g. \code{list(0.4,0.5,0.6)})
#' @param ncol integer, number of columns to have on the plot if image is 3D.
#' @param slice_buffer integer,
#' how many slices to buffer when finding the non-zero slices of a 3D images. So, if \code{slice_buffer=10}, then the first slice in a 3D image will be the first non-zero slice index plus 10 more slices.
#' @param black_bg boolean:
#' if TRUE, the background of the image(s) will be black. if FALSE, the background of the image(s) will be determined by the values `bg_thresh_quant` and `bg_val_quant`.
#' @param bg_thresh_quant float
#' if white_bg=TRUE, the background will be determined by thresholding the image at the `bg_thresh` quantile value and setting the background intensity to the `bg_val` quantile value. This value should be in [0, 1] - somewhere around 0.01 is recommended. - equal to 1 will threshold the entire image - equal to 0 will threshold none of the image
#' @param bg_val_quant float:
#' if white_bg=TRUE, the background will be determined by thresholding the image at the `bg_thresh` quantile value and setting the background intensity to the `bg_val` quantile value. This value should be in [0, 1] - equal to 1 is pure white - equal to 0 is pure black - somewhere in between is gray
#' @param domain_image_map 'ANTsImage':
#' this input 'ANTsImage' or list of 'ANTsImage' types contains a reference image `domain_image` and optional reference mapping named `domainMap`. If supplied, the image(s) to be plotted will be mapped to the domain image space before plotting - useful for non-standard image orientations.
#' @param crop boolean:
#' if true, the image(s) will be cropped to their bounding boxes, resulting in a potentially smaller image size. if false, the image(s) will not be cropped
#' @param scale boolean or 2-tuple:
#' if true, nothing will happen to intensities of image(s) and overlay(s) if false, dynamic range will be maximized when visualizing overlays if 2-list, the image will be dynamically scaled between these quantiles
#' @param reverse boolean:
#' if true, the order in which the slices are plotted will be reversed. This is useful if you want to plot from the front of the brain first to the back of the brain, or vice-versa
#' @param title string:
#' add a title to the plot
#' @param title_fontsize title_fontsize
#' @param title_dx title_dx
#' @param title_dy title_dy
#' @param filename string:
#' if given, the resulting image will be saved to this file
#' @param dpi integer:
#' determines resolution of image if saved to file. Higher values result in higher resolution images, but at a cost of having a larger file size
#' @param figsize figsize
#' @param reorient reorient
#' @param resample bool:
#' if true, resample image if spacing is very unbalanced.
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'   np <- import("numpy")
#'   img <- ants.image_read(ants.get_data('r16'))
#'   segs <- img$kmeans_segmentation(k = 3L)[['segmentation']]
#'
#'   # R cannot call python's `==` or "*", hence use `__ne__` and `__mul__`
#'   overlay <- segs$`__mul__`( segs$`__ne__`(1) )
#'   ants.plot(img, overlay, crop=TRUE)
#' }
#'
#' @export
ants.plot <- function(image, overlay = NULL, blend = FALSE, alpha = 1L, cmap = "Greys_r", overlay_cmap = "turbo", overlay_alpha = 0.9, vminol = NULL, vmaxol = NULL, cbar = FALSE, cbar_length = 0.8, cbar_dx = 0.0, cbar_vertical = TRUE, axis = 0L, nslices = 12L, slices = NULL, ncol = NULL, slice_buffer = NULL, black_bg = TRUE, bg_thresh_quant = 0.01, bg_val_quant = 0.99, domain_image_map = NULL, crop = FALSE, scale = FALSE, reverse = FALSE, title = NULL, title_fontsize = 20L, title_dx = 0.0, title_dy = 0.0, filename = NULL, dpi = 500L, figsize = 1.5, reorient = TRUE, resample = TRUE) {
  ants <- load_ants()
  ants$plot(
    image = image,
    overlay = overlay,
    blend = blend,
    alpha = alpha,
    cmap = cmap,
    overlay_cmap = overlay_cmap,
    overlay_alpha = overlay_alpha,
    vminol = vminol,
    vmaxol = vmaxol,
    cbar = cbar,
    cbar_length = cbar_length,
    cbar_dx = cbar_dx,
    cbar_vertical = cbar_vertical,
    axis = axis,
    nslices = nslices,
    slices = slices,
    ncol = ncol,
    slice_buffer = slice_buffer,
    black_bg = black_bg,
    bg_thresh_quant = bg_thresh_quant,
    bg_val_quant = bg_val_quant,
    domain_image_map = domain_image_map,
    crop = crop,
    scale = scale,
    reverse = reverse,
    title = title,
    title_fontsize = title_fontsize,
    title_dx = title_dx,
    title_dy = title_dy,
    filename = filename,
    dpi = dpi,
    figsize = figsize,
    reorient = reorient,
    resample = resample
  )
}




#' @name ants.plot_directory
#'
#' @title Create and save an 'ANTsPy' plot in a directory
#' @description Create and save an 'ANTsPy' plot for every image matching a
#' given regular expression in a directory optionally recursively. This is a
#' good function for quick visualize exploration of all of images in a directory
#'
#' @param directory string:
#' directory in which to search for images and plot them
#'
#' @param recursive boolean:
#' If true, this function will search through all directories under the given directory recursively to make plots. If false, this function will only create plots for images in the given directory
#'
#' @param regex string:
#' regular expression used to filter out certain filenames or suffixes
#'
#' @param save_prefix string:
#' sub-string that will be appended to the beginning of all saved plot filenames. Default is to add nothing.
#'
#' @param save_suffix string:
#' sub-string that will be appended to the end of all saved plot filenames. Default is add nothing.
#'
#' @param axis whether to plot axis
#' @param ... additional arguments to pass onto the `ants.plot` function. e.g. overlay, alpha, cmap, etc. See `ants.plot` for more options.
#'
#' @examples
#' if(interactive() && ants_available() &&
#'   dir.exists('~/desktop/testdir')) {
#'
#'   ants.plot_directory(directory='~/desktop/testdir',
#'                       recursive=FALSE, regex='*')
#' }
#'
#' @export
ants.plot_directory <- function(directory, recursive = FALSE, regex = "*", save_prefix = "", save_suffix = "", axis = NULL, ...) {
  ants <- load_ants()
  ants$plot_directory(
    directory = directory,
    recursive = recursive,
    regex = regex,
    save_prefix = save_prefix,
    save_suffix = save_suffix,
    axis = axis, ...
  )
}


#' @name ants.plot_grid
#'
#' @title Plot a collection of images in an arbitrarily-defined grid
#'
#' @param figsize,rpad,cpad,vmin,vmax,colorbar,cmap,title,tfontsize,title_dx,title_dy,rlabels,rfontsize,rfontcolor,rfacecolor,clabels,cfontsize,cfontcolor,cfacecolor,filename,dpi,transparent parameters passing to other methods
#'
#' @param images : list of ANTsImage types:
#' image(s) to plot. if one image, this image will be used for all grid locations. if multiple images, they should be arrange in a list the same shape as the `gridsize` argument.
#'
#' @param slices : integer or list of integers:
#' slice indices to plot if one integer, this slice index will be used for all images if multiple integers, they should be arranged in a list the same shape as the `gridsize` argument
#'
#' @param axes : integer or list of integers:
#' axis or axes along which to plot image slices if one integer, this axis will be used for all images if multiple integers, they should be arranged in a list the same shape as the `gridsize` argument
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   np <- import("numpy")
#'
#'   mni1 <- ants.image_read(ants.get_data('mni'))
#'   mni2 <- mni1$smooth_image(1.)
#'   mni3 <- mni1$smooth_image(2.)
#'   mni4 <- mni1$smooth_image(3.)
#'
#'   images <- list(list(mni1, mni2), list(mni3, mni4))
#'   slices <- list(list(100L, 100L), list(100L, 100L))
#'   axes <- list(list(2L, 2L), list(2L, 2L))
#'
#'   # standard plotting
#'   ants.plot_grid(images=images, slices=slices, title='2x2 Grid')
#'
#'   ants.plot_grid(images, slices, cpad=0.02, title='Col Padding')
#'   ants.plot_grid(images, slices, rpad=0.02, title='Row Padding')
#'   ants.plot_grid(images, slices, rpad=0.02, cpad=0.02,
#'                  title='Row and Col Padding')
#'   # Adding plain row and/or column labels
#'   ants.plot_grid(images, slices,
#'                  rlabels=list('Row #1', 'Row #2'))
#'   ants.plot_grid(images, slices,
#'                  clabels=list('Col #1', 'Col #2'))
#'   ants.plot_grid(images, slices,
#'                  rlabels=list('Row 1', 'Row 2'),
#'                  clabels=list('Col 1', 'Col 2'))
#'
#'   images <- list(list(mni1, mni2, mni2),
#'                  list(mni3, mni4, mni4))
#'   slices <- list(list(100L, 100L, 100L),
#'                  list(100L, 100L, 100L))
#'   axes <- list(list(0L, 1L, 2L),
#'                list(0L, 1L, 2L))
#'   ants.plot_grid(images, slices, axes,
#'                  title='Publication Figures with ANTsPy',
#'                  tfontsize=20, title_dy=0.03, title_dx=-0.04,
#'                  rlabels=list('Row 1', 'Row 2'),
#'                  clabels=list('Col 1', 'Col 2', 'Col 3'),
#'                  rfontsize=16, cfontsize=16)
#' }
#'
#'
#' @export
ants.plot_grid <- function(images, slices = NULL, axes = 2L, figsize = 1.0, rpad = 0L, cpad = 0L, vmin = NULL, vmax = NULL, colorbar = TRUE, cmap = "Greys_r", title = NULL, tfontsize = 20L, title_dx = 0L, title_dy = 0L, rlabels = NULL, rfontsize = 14L, rfontcolor = "white", rfacecolor = "black", clabels = NULL, cfontsize = 14L, cfontcolor = "white", cfacecolor = "black", filename = NULL, dpi = 400L, transparent = TRUE, ...) {
  ants <- load_ants()
  ants$plot_grid(
    images = images,
    slices = slices,
    axes = axes,
    figsize = figsize,
    rpad = rpad,
    cpad = cpad,
    vmin = vmin,
    vmax = vmax,
    colorbar = colorbar,
    cmap = cmap,
    title = title,
    tfontsize = tfontsize,
    title_dx = title_dx,
    title_dy = title_dy,
    rlabels = rlabels,
    rfontsize = rfontsize,
    rfontcolor = rfontcolor,
    rfacecolor = rfacecolor,
    clabels = clabels,
    cfontsize = cfontsize,
    cfontcolor = cfontcolor,
    cfacecolor = cfacecolor,
    filename = filename,
    dpi = dpi,
    transparent = transparent,
    ...
  )
}


#' @name ants.plot_hist
#' @title Plot a histogram from an ANTsImage
#' @param image ANTsImage image from which histogram will be created
#' @param threshold threshold
#' @param fit_line fit_line
#' @param normfreq normfreq
#' @param title title
#' @param grid grid
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param facecolor facecolor
#' @param alpha alpha
#'
#' @export
ants.plot_hist <- function(image, threshold = 0.0, fit_line = FALSE, normfreq = TRUE, title = NULL, grid = TRUE, xlabel = NULL, ylabel = NULL, facecolor = "green", alpha = 0.75) {
  ants <- load_ants()
  ants$plot_hist(
    image = image,
    threshold = threshold,
    fit_line = fit_line,
    normfreq = normfreq,
    title = title,
    grid = grid,
    xlabel = xlabel,
    ylabel = ylabel,
    facecolor = facecolor,
    alpha = alpha
  )
}


#' @name ants.plot_ortho
#' @title Plot an orthographic view of a 3D image
#' @description
#' Use mask_image and/or threshold_image to preprocess images to be be overlaid and display the overlays in a given range.
#' See the wiki examples. TODO: - add colorbar option ANTsR function: N/A Arguments ---------
#'
#' @param image : ANTsImage image to plot
#' @param overlay : ANTsImage image to overlay on base image
#' @param slices : list or list of 3 integers slice indices along each axis to plot. This can be absolute array indices (e.g. (80,100,120)), or this can be relative array indices (e.g. (0.4,0.5,0.6)). The default is to take the middle slice along each axis.
#' @param xyz : list or list of 3 integers if given, solid lines will be drawn to converge at this coordinate. This is useful for pinpointing a specific location in the image.
#' @param flat : boolean if true, the ortho image will be plot in one row if false, the ortho image will be a 2x2 grid with the bottom left corner blank
#' @param cmap : string colormap to use for base image. See matplotlib.
#' @param overlay_cmap : string colormap to use for overlay images, if applicable. See matplotlib.
#' @param overlay_alpha : float level of transparency for any overlays. Smaller value means the overlay is more transparent. See matplotlib.
#' @param axis : integer which axis to plot along if image is 3D
#' @param black_bg : boolean if TRUE, the background of the image(s) will be black. if FALSE, the background of the image(s) will be determined by the values `bg_thresh_quant` and `bg_val_quant`.
#' @param bg_thresh_quant : float if white_bg=TRUE, the background will be determined by thresholding the image at the `bg_thresh` quantile value and setting the background intensity to the `bg_val` quantile value. This value should be in [0, 1] - somewhere around 0.01 is recommended. - equal to 1 will threshold the entire image - equal to 0 will threshold none of the image
#' @param bg_val_quant : float if white_bg=TRUE, the background will be determined by thresholding the image at the `bg_thresh` quantile value and setting the background intensity to the `bg_val` quantile value. This value should be in [0, 1] - equal to 1 is pure white - equal to 0 is pure black - somewhere in between is gray
#' @param domain_image_map : ANTsImage this input ANTsImage or list of ANTsImage types contains a reference image `domain_image` and optional reference mapping named `domainMap`. If supplied, the image(s) to be plotted will be mapped to the domain image space before plotting - useful for non-standard image orientations.
#' @param crop : boolean if true, the image(s) will be cropped to their bounding boxes, resulting in a potentially smaller image size. if false, the image(s) will not be cropped
#' @param scale : boolean or 2-list if true, nothing will happen to intensities of image(s) and overlay(s) if false, dynamic range will be maximized when visualizing overlays if 2-list, the image will be dynamically scaled between these quantiles
#' @param title : string add a title to the plot
#' @param filename : string if given, the resulting image will be saved to this file
#' @param dpi : integer determines resolution of image if saved to file. Higher values result in higher resolution images, but at a cost of having a larger file size
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   mni = ants.image_read(ants.get_data('mni'))
#'   ants.plot_ortho(mni, xyz=list(100L,100L,100L))
#'
#'   mni2 = mni$threshold_image(7000, mni$max())
#'   ants.plot_ortho(mni, overlay=mni2)
#'   ants.plot_ortho(mni, overlay=mni2, flat=TRUE)
#'
#'   ants.plot_ortho(mni, overlay=mni2, xyz=list(110L,110L,110L),
#'                   xyz_lines=FALSE, text='Lines Turned Off', textfontsize=12L)
#'
#'   ants.plot_ortho(mni, mni2, xyz=list(120L,100L,100L),
#'                   text=' Example Ortho Text', textfontsize=12L,
#'                   title='Example Ortho Title', titlefontsize=12L)
#'
#' }
#'
#' @export
ants.plot_ortho <- function(image, overlay = NULL, reorient = TRUE, blend = FALSE, xyz = NULL, xyz_lines = TRUE, xyz_color = "red", xyz_alpha = 0.6, xyz_linewidth = 2L, xyz_pad = 5L, orient_labels = TRUE, alpha = 1L, cmap = "Greys_r", overlay_cmap = "jet", overlay_alpha = 0.9, black_bg = TRUE, bg_thresh_quant = 0.01, bg_val_quant = 0.99, crop = FALSE, scale = FALSE, domain_image_map = NULL, title = NULL, titlefontsize = 24L, title_dx = 0L, title_dy = 0L, text = NULL, textfontsize = 24L, textfontcolor = "white", text_dx = 0L, text_dy = 0L, filename = NULL, dpi = 500L, figsize = 1.0, flat = FALSE, transparent = TRUE) {
  ants <- load_ants()
  ants$plot_ortho(
    image = image,
    overlay = overlay,
    reorient = reorient,
    blend = blend,
    xyz = xyz,
    xyz_lines = xyz_lines,
    xyz_color = xyz_color,
    xyz_alpha = xyz_alpha,
    xyz_linewidth = xyz_linewidth,
    xyz_pad = xyz_pad,
    orient_labels = orient_labels,
    alpha = alpha,
    cmap = cmap,
    overlay_cmap = overlay_cmap,
    overlay_alpha = overlay_alpha,
    black_bg = black_bg,
    bg_thresh_quant = bg_thresh_quant,
    bg_val_quant = bg_val_quant,
    crop = crop,
    scale = scale,
    domain_image_map = domain_image_map,
    title = title,
    titlefontsize = titlefontsize,
    title_dx = title_dx,
    title_dy = title_dy,
    text = text,
    textfontsize = textfontsize,
    textfontcolor = textfontcolor,
    text_dx = text_dx,
    text_dy = text_dy,
    filename = filename,
    dpi = dpi,
    figsize = figsize,
    flat = flat,
    transparent = transparent
  )
}


#' @rdname ants.plot_ortho
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   mni <- ants.image_read(ants.get_data('mni'))
#'   ch2 <- ants.image_read(ants.get_data('ch2'))
#'   ants.plot_ortho_double(mni, ch2)
#'
#' }
#'
#' @export
ants.plot_ortho_double <- function(image, image2, overlay = NULL, overlay2 = NULL, reorient = TRUE, xyz = NULL, xyz_lines = TRUE, xyz_color = "red", xyz_alpha = 0.6, xyz_linewidth = 2L, xyz_pad = 5L, cmap = "Greys_r", alpha = 1L, cmap2 = "Greys_r", alpha2 = 1L, overlay_cmap = "jet", overlay_alpha = 0.9, overlay_cmap2 = "jet", overlay_alpha2 = 0.9, black_bg = TRUE, bg_thresh_quant = 0.01, bg_val_quant = 0.99, crop = FALSE, scale = FALSE, crop2 = FALSE, scale2 = TRUE, domain_image_map = NULL, title = NULL, titlefontsize = 24L, title_dx = 0L, title_dy = 0L, text = NULL, textfontsize = 24L, textfontcolor = "white", text_dx = 0L, text_dy = 0L, filename = NULL, dpi = 500L, figsize = 1.0, flat = TRUE, transpose = FALSE, transparent = TRUE) {
  ants <- load_ants()
  ants$plot_ortho_double(
    image = image,
    image2 = image2,
    overlay = overlay,
    overlay2 = overlay2,
    reorient = reorient,
    xyz = xyz,
    xyz_lines = xyz_lines,
    xyz_color = xyz_color,
    xyz_alpha = xyz_alpha,
    xyz_linewidth = xyz_linewidth,
    xyz_pad = xyz_pad,
    cmap = cmap,
    alpha = alpha,
    cmap2 = cmap2,
    alpha2 = alpha2,
    overlay_cmap = overlay_cmap,
    overlay_alpha = overlay_alpha,
    overlay_cmap2 = overlay_cmap2,
    overlay_alpha2 = overlay_alpha2,
    black_bg = black_bg,
    bg_thresh_quant = bg_thresh_quant,
    bg_val_quant = bg_val_quant,
    crop = crop,
    scale = scale,
    crop2 = crop2,
    scale2 = scale2,
    domain_image_map = domain_image_map,
    title = title,
    titlefontsize = titlefontsize,
    title_dx = title_dx,
    title_dy = title_dy,
    text = text,
    textfontsize = textfontsize,
    textfontcolor = textfontcolor,
    text_dx = text_dx,
    text_dy = text_dy,
    filename = filename,
    dpi = dpi,
    figsize = figsize,
    flat = flat,
    transpose = transpose,
    transparent = transparent
  )
}


#' @rdname ants.plot_ortho
#'
#' @examples
#'
#' if(interactive() && ants_available()) {
#'
#'   mni = ants.image_read(ants.get_data('mni'))
#'   ch2 = ants.image_read(ants.get_data('ch2'))
#'
#'   ants.plot_ortho_stack( list(mni,mni,mni) )
#'
#' }
#'
#' @export
ants.plot_ortho_stack <- function(images, overlays = NULL, reorient = TRUE, xyz = NULL, xyz_lines = FALSE, xyz_color = "red", xyz_alpha = 0.6, xyz_linewidth = 2L, xyz_pad = 5L, cmap = "Greys_r", alpha = 1L, overlay_cmap = "jet", overlay_alpha = 0.9, black_bg = TRUE, bg_thresh_quant = 0.01, bg_val_quant = 0.99, crop = FALSE, scale = FALSE, domain_image_map = NULL, title = NULL, titlefontsize = 24L, title_dx = 0L, title_dy = 0L, text = NULL, textfontsize = 24L, textfontcolor = "white", text_dx = 0L, text_dy = 0L, filename = NULL, dpi = 500L, figsize = 1.0, colpad = 0L, rowpad = 0L, transpose = FALSE, transparent = TRUE, orient_labels = TRUE) {
  ants <- load_ants()
  ants$plot_ortho_stack(
    images = images,
    overlays = overlays,
    reorient = reorient,
    xyz = xyz,
    xyz_lines = xyz_lines,
    xyz_color = xyz_color,
    xyz_alpha = xyz_alpha,
    xyz_linewidth = xyz_linewidth,
    xyz_pad = xyz_pad,
    cmap = cmap,
    alpha = alpha,
    overlay_cmap = overlay_cmap,
    overlay_alpha = overlay_alpha,
    black_bg = black_bg,
    bg_thresh_quant = bg_thresh_quant,
    bg_val_quant = bg_val_quant,
    crop = crop,
    scale = scale,
    domain_image_map = domain_image_map,
    title = title,
    titlefontsize = titlefontsize,
    title_dx = title_dx,
    title_dy = title_dy,
    text = text,
    textfontsize = textfontsize,
    textfontcolor = textfontcolor,
    text_dx = text_dx,
    text_dy = text_dy,
    filename = filename,
    dpi = dpi,
    figsize = figsize,
    colpad = colpad,
    rowpad = rowpad,
    transpose = transpose,
    transparent = transparent,
    orient_labels = orient_labels
  )
}
