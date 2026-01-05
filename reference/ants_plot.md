# Plot single `'ANTsImage'`

Plot single `'ANTsImage'`

## Usage

``` r
ants_plot(
  image,
  overlay = NULL,
  blend = FALSE,
  alpha = 1,
  cmap = "Greys_r",
  overlay_cmap = "turbo",
  overlay_alpha = 0.9,
  vminol = NULL,
  vmaxol = NULL,
  cbar = FALSE,
  cbar_length = 0.8,
  cbar_dx = 0,
  cbar_vertical = TRUE,
  axis = 0,
  nslices = 12,
  slices = NULL,
  ncol = NULL,
  slice_buffer = NULL,
  black_bg = TRUE,
  bg_thresh_quant = 0.01,
  bg_val_quant = 0.99,
  domain_image_map = NULL,
  crop = FALSE,
  scale = FALSE,
  reverse = FALSE,
  title = "",
  title_fontsize = 20,
  title_dx = 0,
  title_dy = 0,
  filename = NULL,
  dpi = 500,
  figsize = 1.5,
  reorient = TRUE,
  resample = TRUE,
  force_agg = FALSE,
  close_figure = TRUE
)
```

## Arguments

- image:

  `'ANTsImage'`, or something can be converted to `'ANTsImage'`

- overlay:

  overlay `'ANTsImage'`, can be `NULL`, optional

- blend:

  whether to blend image with overlay; default is false

- cmap, alpha:

  image color map and transparency

- overlay_cmap, overlay_alpha:

  overlay color map and transparency

- vminol, vmaxol:

  I could not find its usage

- cbar:

  whether to draw color legend

- cbar_length, cbar_dx, cbar_vertical:

  legend position and size

- axis:

  see 'Details'

- nslices, slices, ncol:

  controls slice to show

- slice_buffer:

  performance

- black_bg, bg_thresh_quant, bg_val_quant:

  controls background

- domain_image_map:

  optional `'ANTsImage'`

- crop, scale, reverse:

  whether to crop, scale, or reverse the image according to background

- title, title_fontsize, title_dx, title_dy:

  image title

- filename, dpi, figsize:

  needed when saving to file

- reorient:

  whether to reorient to `'LAI'` before plotting; default is true

- resample:

  whether to resample

- force_agg:

  whether to force graphic engine to use `'agg'` device; default is
  false

- close_figure:

  whether to close figure when returning the function

## Value

Nothing

## Details

By default, images will be reoriented to `'LAI'` orientation before
plotting. So, if `axis=0`, the images will be ordered from the left side
of the brain to the right side of the brain. If `axis=1`, the images
will be ordered from the anterior (front) of the brain to the posterior
(back) of the brain. And if `axis=2`, the images will be ordered from
the inferior (bottom) of the brain to the superior (top) of the brain.

## Examples

``` r

if(interactive() && ants_available()) {
  ants <- load_ants()
  img <- ants$image_read(ants$get_ants_data('mni'))

  ants_plot(
    img, nslices = 12, black_bg = FALSE,
    bg_thresh_quant = 0.05, bg_val_quant = 1.0, axis = 2,
    cbar = TRUE, crop = TRUE, reverse = TRUE, cbar_vertical = FALSE,
    ncol = 4, title = "Axial view of MNI brain"
  )
}

```
