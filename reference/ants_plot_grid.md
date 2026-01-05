# Plot multiple `'ANTsImage'`

R-friendly wrapper function for `ants$plot_grid`

## Usage

``` r
ants_plot_grid(
  images,
  shape = NULL,
  slices = 0,
  axes = 2,
  figsize = 1,
  rpad = 0,
  cpad = 0,
  vmin = NULL,
  vmax = NULL,
  colorbar = TRUE,
  cmap = "Greys_r",
  title = "",
  tfontsize = 20,
  title_dx = 0,
  title_dy = 0,
  rlabels = NULL,
  rfontsize = 14,
  rfontcolor = "black",
  rfacecolor = "white",
  clabels = NULL,
  cfontsize = 14,
  cfontcolor = "black",
  cfacecolor = "white",
  filename = NULL,
  dpi = 400,
  transparent = TRUE,
  ...,
  force_agg = FALSE,
  close_figure = TRUE
)
```

## Arguments

- images:

  a single `'ANTsImage'`, list, or nested list of `'ANTsImage'`

- shape:

  shape of grid, default is using dimensions of `images`

- slices:

  length of one or equaling to length of `slices`, slice number to plot

- axes:

  `0` for `'sagittal'`, `1` for `'coronal'`, `2` for `'axial'`; default
  is `2`

- figsize, rpad, cpad, colorbar, cmap, transparent:

  graphical parameters

- vmin, vmax:

  value threshold for the image

- title:

  title of figure

- title_dx, title_dy, tfontsize:

  controls title margin and size

- rlabels, clabels:

  row and column labels

- rfontsize, rfontcolor, rfacecolor, cfontsize, cfontcolor, cfacecolor:

  row and column font size, color, and background color

- filename, dpi:

  parameters to save figures

- ...:

  passed to `ants$plot_grid`; make sure all entries are named

- force_agg:

  whether to force graphic engine to use `'agg'` device; default is
  false

- close_figure:

  whether to close figure when returning the function

## Value

Nothing

## Examples

``` r
if(interactive() && ants_available()) {
  ants <- load_ants()
  image1 <- ants$image_read(ants$get_ants_data('mni'))
  image2 <- image1$smooth_image(1.0)
  image3 <- image1$smooth_image(2.0)
  image4 <- image1$smooth_image(3.0)

  ants_plot_grid(
    list(image1, image2, image3, image4),
    slices = 100, title = "4x1 Grid"
  )

  ants_plot_grid(
    list(image1, image2, image3, image4),
    shape = c(2, 2),
    slices = 100, title = "2x2 Grid"
  )
  ants_plot_grid(
    list(image1, image2, image3, image4),
    shape = c(2, 2), axes = c(0,1,2,1),
    slices = 100, title = "2x2 Grid (diff. anatomical slices)"
  )

}


```
