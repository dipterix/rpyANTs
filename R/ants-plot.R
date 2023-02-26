#' @title Plot single \code{'ANTsImage'}
#' @param image \code{'ANTsImage'}, or something can be converted to \code{'ANTsImage'}
#' @param overlay overlay \code{'ANTsImage'}, can be \code{NULL}, optional
#' @param blend whether to blend image with overlay; default is false
#' @param cmap,alpha image color map and transparency
#' @param overlay_cmap,overlay_alpha overlay color map and transparency
#' @param vminol,vmaxol I could not find its usage
#' @param cbar whether to draw color legend
#' @param cbar_length,cbar_dx,cbar_vertical legend position and size
#' @param axis see 'Details'
#' @param nslices,slices,ncol controls slice to show
#' @param slice_buffer performance
#' @param black_bg,bg_thresh_quant,bg_val_quant controls background
#' @param domain_image_map optional \code{'ANTsImage'}
#' @param crop,scale,reverse whether to crop, scale, or reverse the image
#' according to background
#' @param title,title_fontsize,title_dx,title_dy image title
#' @param filename,dpi,figsize needed when saving to file
#' @param reorient whether to reorient to \code{'LAI'} before plotting;
#' default is true
#' @param resample whether to resample
#' @param force_agg whether to force graphic engine to use \code{'agg'} device;
#' default is false
#' @param close_figure whether to close figure when returning the function
#' @returns Nothing
#'
#' @details
#' By default, images will be reoriented to \code{'LAI'} orientation before
#' plotting. So, if \code{axis=0}, the images will be ordered from the
#' left side of the brain to the right side of the brain. If \code{axis=1},
#' the images will be ordered from the anterior (front) of the brain to
#' the posterior (back) of the brain. And if \code{axis=2}, the images will
#' be ordered from the inferior (bottom) of the brain to the superior (top)
#' of the brain.
#'
#' @examples
#'
#'
#' if(interactive() && ants_available()) {
#'   ants <- load_ants()
#'   img <- ants$image_read(ants$get_ants_data('mni'))
#'
#'   ants_plot(
#'     img, nslices = 12, black_bg = FALSE,
#'     bg_thresh_quant = 0.05, bg_val_quant = 1.0, axis = 2,
#'     cbar = TRUE, crop = TRUE, reverse = TRUE, cbar_vertical = FALSE,
#'     ncol = 4, title = "Axial view of MNI brain"
#'   )
#' }
#'
#'
#' @export
ants_plot <- function(
    image, overlay=NULL, blend=FALSE, alpha=1,
    cmap='Greys_r', overlay_cmap='turbo',
    overlay_alpha=0.9, vminol=NULL, vmaxol=NULL,
    cbar=FALSE, cbar_length=0.8, cbar_dx=0.0, cbar_vertical=TRUE,
    axis=0, nslices=12, slices=NULL, ncol=NULL, slice_buffer=NULL,
    black_bg=TRUE, bg_thresh_quant=0.01, bg_val_quant=0.99,
    domain_image_map=NULL, crop=FALSE, scale=FALSE, reverse=FALSE,
    title="", title_fontsize=20, title_dx=0.0, title_dy=0.0,
    filename=NULL, dpi=500, figsize=1.5, reorient=TRUE, resample=TRUE,
    force_agg = FALSE, close_figure = TRUE
  ) {

  # DIPSAUS DEBUG START
  # image <- '~/.antspy/r16slice.jpg'
  # list2env(list(overlay=NULL, blend=FALSE, alpha=1,
  #               cmap='Greys_r', overlay_cmap='turbo',
  #               overlay_alpha=0.4, vminol=NULL, vmaxol=NULL,
  #               cbar=FALSE, cbar_length=0.8, cbar_dx=0.0, cbar_vertical=TRUE,
  #               axis=0, nslices=12, slices=NULL, ncol=NULL, slice_buffer=NULL,
  #               black_bg=TRUE, bg_thresh_quant=0.01, bg_val_quant=0.99,
  #               domain_image_map=NULL, crop=FALSE, scale=FALSE, reverse=FALSE,
  #               title="", title_fontsize=20, title_dx=0.0, title_dy=0.0,
  #               filename=NULL, dpi=500, figsize=1.5, reorient=TRUE, resample=TRUE), envir=.GlobalEnv)


  blend <- convert_if_not_python(blend, { as.logical(blend) })
  cbar <- convert_if_not_python(cbar, { as.logical(cbar) })
  cbar_vertical <- convert_if_not_python(cbar_vertical, { as.logical(cbar_vertical) })
  black_bg <- convert_if_not_python(black_bg, { as.logical(black_bg) })
  crop <- convert_if_not_python(crop, { as.logical(crop) })
  scale <- convert_if_not_python(scale, { as.logical(scale) })
  reverse <- convert_if_not_python(reverse, { as.logical(reverse) })
  reorient <- convert_if_not_python(reorient, { as.logical(reorient) })
  resample <- convert_if_not_python(resample, { as.logical(resample) })

  alpha <- convert_if_not_python(alpha, { as.double(alpha) })
  overlay_alpha <- convert_if_not_python(overlay_alpha, { as.double(overlay_alpha) })
  title_dx <- convert_if_not_python(title_dx, { as.double(title_dx) })
  title_dy <- convert_if_not_python(title_dy, { as.double(title_dy) })
  bg_thresh_quant <- convert_if_not_python(bg_thresh_quant, { as.double(bg_thresh_quant) })
  bg_val_quant <- convert_if_not_python(bg_val_quant, { as.double(bg_val_quant) })
  figsize <- convert_if_not_python(figsize, { as.double(figsize) })
  cbar_length <- convert_if_not_python(cbar_length, { as.double(cbar_length) })
  cbar_dx <- convert_if_not_python(cbar_dx, { as.double(cbar_dx) })
  vminol <- convert_if_not_python(vminol, {
    if(length(vminol)) { as.double(vminol) } else {NULL }
  })
  vmaxol <- convert_if_not_python(vmaxol, {
    if(length(vmaxol)) { as.double(vmaxol) } else {NULL }
  })

  cmap <- convert_if_not_python(cmap, { as.character(cmap) })
  overlay_cmap <- convert_if_not_python(overlay_cmap, { as.character(overlay_cmap) })

  axis <- convert_if_not_python(axis, { as.integer(axis) })
  nslices <- convert_if_not_python(nslices, { as.integer(nslices) })
  ncol <- convert_if_not_python(ncol, { as.integer(ncol) })
  title_fontsize <- convert_if_not_python(title_fontsize, { as.integer(title_fontsize) })
  dpi <- convert_if_not_python(dpi, { as.integer(dpi) })
  slice_buffer <- convert_if_not_python(slice_buffer, {
    if(length(slice_buffer)) { as.integer(slice_buffer) } else { NULL }
  })

  slices <- convert_if_not_python(slices, {
    if(length(slices)) {
      slices <- as.numeric(slices)
      if(any(slices < 1 & slices > 0)) {
        slices <- as.double(slices)
      } else {
        slices <- as.integer(slices)
      }
      slices
    } else {
      NULL
    }
  })
  title <- convert_if_not_python(title, { paste(title, collapse = "") })
  filename <- convert_if_not_python(filename, {
    if(length(filename) != 1) {
      NULL
    } else {
      filename
    }
  })

  image <- as_ANTsImage(image, strict = TRUE)
  overlay <- as_ANTsImage(overlay, strict = FALSE)
  domain_image_map <- as_ANTsImage(domain_image_map, strict = FALSE)


  ants <- load_ants()
  matplotlib <- import("matplotlib", convert = FALSE)
  if( force_agg ) {
    matplotlib$use("Agg", force = TRUE)
  }

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
    ncol = ncol,
    cbar_length = cbar_length,
    cbar_dx = cbar_dx,
    cbar_vertical = cbar_vertical,
    axis = axis,
    nslices = nslices,
    slices = slices,
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
    resample= resample
  )


  if(close_figure) {
    matplotlib$pyplot$close()
  }
  return(invisible())
}



#' @title Plot multiple \code{'ANTsImage'}
#' @description
#' R-friendly wrapper function for \code{ants$plot_grid}
#' @param images a single \code{'ANTsImage'}, list, or nested list of \code{'ANTsImage'}
#' @param shape shape of grid, default is using dimensions of \code{images}
#' @param slices length of one or equaling to length of \code{slices}, slice
#' number to plot
#' @param axes \code{0} for \code{'sagittal'}, \code{1} for \code{'coronal'},
#' \code{2} for \code{'axial'}; default is \code{2}
#' @param vmin,vmax value threshold for the image
#' @param title title of figure
#' @param title_dx,title_dy,tfontsize controls title margin and size
#' @param figsize,rpad,cpad,colorbar,cmap,transparent graphical parameters
#' @param rlabels,clabels row and column labels
#' @param rfontsize,rfontcolor,rfacecolor,cfontsize,cfontcolor,cfacecolor row
#' and column font size, color, and background color
#' @param filename,dpi parameters to save figures
#' @param force_agg whether to force graphic engine to use \code{'agg'} device;
#' default is false
#' @param close_figure whether to close figure when returning the function
#' @param ... passed to \code{ants$plot_grid}; make sure all entries are named
#' @returns Nothing
#' @examples
#'
#' if(interactive() && ants_available()) {
#'   ants <- load_ants()
#'   image1 <- ants$image_read(ants$get_ants_data('mni'))
#'   image2 <- image1$smooth_image(1.0)
#'   image3 <- image1$smooth_image(2.0)
#'   image4 <- image1$smooth_image(3.0)
#'
#'   ants_plot_grid(
#'     list(image1, image2, image3, image4),
#'     slices = 100, title = "4x1 Grid"
#'   )
#'
#'   ants_plot_grid(
#'     list(image1, image2, image3, image4),
#'     shape = c(2, 2),
#'     slices = 100, title = "2x2 Grid"
#'   )
#'   ants_plot_grid(
#'     list(image1, image2, image3, image4),
#'     shape = c(2, 2), axes = c(0,1,2,1),
#'     slices = 100, title = "2x2 Grid (diff. anatomical slices)"
#'   )
#'
#' }
#'
#'
#'
#' @export
ants_plot_grid <- function(
    images, shape = NULL, slices=0, axes=2, figsize=1.0, rpad=0, cpad=0,
    vmin=NULL, vmax=NULL, colorbar=TRUE, cmap='Greys_r', title="",
    tfontsize=20, title_dx=0, title_dy=0, rlabels=NULL, rfontsize=14,
    rfontcolor='black', rfacecolor='white', clabels=NULL,
    cfontsize=14, cfontcolor='black', cfacecolor='white',
    filename=NULL, dpi=400, transparent=TRUE, ...,
    force_agg = FALSE, close_figure = TRUE) {

  # DIPSAUS DEBUG START
  # mni1 = ants$image_read(ants$get_data('mni'))
  # images <- list(mni1, mni1, mni1, mni1)
  # list2env(list(axes=2, figsize=1.0, rpad=0, cpad=0,
  #               vmin=NULL, vmax=NULL, colorbar=TRUE, cmap='Greys_r', title="",
  #               tfontsize=20, title_dx=0, title_dy=0, rlabels=NULL, rfontsize=14,
  #               rfontcolor='white', rfacecolor='black', clabels=NULL,
  #               cfontsize=14, cfontcolor='white', cfacecolor='black',
  #               filename=NULL, dpi=400, transparent=TRUE), envir=.GlobalEnv)
  # slices <- 100



  np <- import("numpy", convert = FALSE)
  ants <- load_ants()

  if(inherits(images, "ants.core.ants_image.ANTsImage")) {
    images <- list(images)
  }
  images <- np$asarray(images)

  if(length(shape)) {
    images <- images$reshape(tuple(as.list(as.integer(shape))))
  }
  if(py_len(images$shape) == 1) {
    n <- c(1L, as.integer(to_r(images$shape)))
    images <- images$reshape(n)
  }

  slices <- convert_if_not_python(slices, {
    if(!length(slices)) {
      slices <- 0L
    }
    if(length(slices) == 1) {
      slices <- rep(as.integer(slices), prod(as.integer(to_r(images$shape))))
    } else {
      slices <- np$asarray(as.integer(unlist(slices)), dtype = "int")
    }
  })

  if(!inherits(slices, "python.builtin.NoneType")) {
    if(!inherits(slices, c("python.builtin.int", "python.builtin.float"))) {
      if(!inherits(slices, "numpy.ndarray")) {
        slices <- np$asarray(slices, dtype = "int")
      }
      slices <- slices$reshape(images$shape)
    }
  }

  axes <- convert_if_not_python(axes, as.integer(axes))
  if(length(axes) > 1) {
    axes <- np$asarray(axes)
    axes <- axes$reshape(images$shape)
  }

  figsize <- convert_if_not_python(figsize, as.double(figsize))
  rpad <- convert_if_not_python(rpad, as.double(rpad))
  cpad <- convert_if_not_python(cpad, as.double(cpad))
  vmin <- convert_if_not_python(vmin, if(length(vmin)) { as.double(vmin) } else { py_none() })
  vmax <- convert_if_not_python(vmax, if(length(vmax)) { as.double(vmax) } else { py_none() })
  colorbar <- convert_if_not_python(colorbar, isTRUE(as.logical(colorbar)))
  cmap <- convert_if_not_python(cmap, paste(as.character(cmap), collapse = ""))
  title <- convert_if_not_python(title, paste(as.character(title), collapse = ""))
  tfontsize <- convert_if_not_python(tfontsize, as.double(tfontsize))
  title_dx <- convert_if_not_python(title_dx, as.double(title_dx))
  title_dy <- convert_if_not_python(title_dy, as.double(title_dy))
  rlabels <- convert_if_not_python(rlabels, {
    if(length(rlabels)) { py_list(as.list(as.character(rlabels))) } else { py_none() }
  })
  clabels <- convert_if_not_python(clabels, {
    if(length(clabels)) { py_list(as.list(as.character(clabels))) } else { py_none() }
  })
  rfontsize <- convert_if_not_python(rfontsize, as.double(rfontsize))
  rfontcolor <- convert_if_not_python(rfontcolor, as_hexcolor(rfontcolor))
  rfacecolor <- convert_if_not_python(rfacecolor, as_hexcolor(rfacecolor))
  cfontsize <- convert_if_not_python(rfacecolor, as.double(cfontsize))
  cfontcolor <- convert_if_not_python(cfontcolor, as_hexcolor(cfontcolor))
  cfacecolor <- convert_if_not_python(cfacecolor, as_hexcolor(cfacecolor))

  filename <- convert_if_not_python(filename, {
    if(length(filename)) {
      paste(as.character(filename), collapse = "")
    } else {
      NULL
    }
  })
  dpi <- convert_if_not_python(dpi, as.integer(dpi))
  transparent <- convert_if_not_python(transparent, isTRUE(as.logical(transparent)))

  matplotlib <- import("matplotlib", convert = FALSE)
  if( force_agg ) {
    matplotlib$use("Agg", force = TRUE)
  }

  ants$plot_grid(images=images, slices=slices, axes = axes, figsize = figsize,
                 rpad = rpad, cpad = cpad, vmin = vmin, vmax = vmax,
                 colorbar = colorbar, cmap = cmap, title = title,
                 tfontsize = tfontsize, title_dx = title_dx,
                 rlabels = rlabels, clabels = clabels,
                 rfontsize = rfontsize, rfontcolor = rfontcolor, rfacecolor = rfacecolor,
                 cfontsize = cfontsize, cfontcolor = cfontcolor, cfacecolor = cfacecolor,
                 filename = filename, dpi = dpi, transparent = transparent)

  if(close_figure) {
    matplotlib$pyplot$close()
  }


  return(invisible())

}
