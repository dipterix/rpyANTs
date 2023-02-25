#' @title Plot multiple 'ANTsImage'
#' @description
#' R-friendly wrapper function for \code{ants$plot_grid}
#' @param images a single 'ANTsImage', list, or nested list of 'ANTsImage'
#' @param shape shape of grid, default is using dimensions of \code{images}
#' @param slices length of one or equaling to length of \code{slices}, slice
#' number to plot
#' @param axes \code{0} for 'sagittal', \code{1} for 'coronal',
#' \code{2} for 'axial'; default is \code{2}
#' @param vmin,vmax value threshold for the image
#' @param title title of figure
#' @param title_dx,title_dy,tfontsize controls title margin and size
#' @param figsize,rpad,cpad,colorbar,cmap,transparent graphical parameters
#' @param rlabels,clabels row and column labels
#' @param rfontsize,rfontcolor,rfacecolor,cfontsize,cfontcolor,cfacecolor row
#' and column font size, color, and background color
#' @param filename,dpi parameters to save figures
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
    filename=NULL, dpi=400, transparent=TRUE, ...) {

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
  matplotlib <- import("matplotlib", convert = FALSE)
  matplotlib$use("Agg", force = TRUE)
  ants <- load_ants()

  if(inherits(images, "ants.core.ants_image.ANTsImage")) {
    images <- list(images)
  }
  images <- np$asarray(images)

  if(length(shape)) {
    images <- images$reshape(tuple(as.list(as.integer(shape))))
  }
  if(py_len(images$shape) == 1) {
    n <- c(1L, as.integer(py_to_r(images$shape)))
    images <- images$reshape(n)
  }

  slices <- convert_if_not_python(slices, {
    if(!length(slices)) {
      slices <- 0L
    }
    if(length(slices) == 1) {
      slices <- rep(as.integer(slices), prod(as.integer(py_to_r(images$shape))))
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

  ants$plot_grid(images=images, slices=slices, axes = axes, figsize = figsize,
                 rpad = rpad, cpad = cpad, vmin = vmin, vmax = vmax,
                 colorbar = colorbar, cmap = cmap, title = title,
                 tfontsize = tfontsize, title_dx = title_dx,
                 rlabels = rlabels, clabels = clabels,
                 rfontsize = rfontsize, rfontcolor = rfontcolor, rfacecolor = rfacecolor,
                 cfontsize = cfontsize, cfontcolor = cfontcolor, cfacecolor = cfacecolor,
                 filename = filename, dpi = dpi, transparent = transparent)
  matplotlib$pyplot$close()

}
