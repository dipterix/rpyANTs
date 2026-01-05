# Process brain image prior to segmentation

Strip skulls, normalize intensity, align and re-sample to template. This
procedure is needed for many `antspynet` functions since the deep neural
networks are trained in template spaces

## Usage

``` r
antspynet_preprocess_brain_image(
  x,
  truncate_intensity = c(0.01, 0.99),
  brain_extraction_modality = c("none", "t1", "t1v0", "t1nobrainer", "t1combined",
    "flair", "t2", "bold", "fa", "t1infant", "t2infant"),
  template_transform_type = c("None", "Affine", "Rigid"),
  template = c("biobank", "croppedMni152"),
  do_bias_correction = TRUE,
  return_bias_field = FALSE,
  do_denoising = TRUE,
  intensity_matching_type = c("regression", "histogram"),
  reference_image = NULL,
  intensity_normalization_type = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  `'ANTsImage'` or path to image to process

- truncate_intensity:

  defines the quantile threshold for truncating the image intensity

- brain_extraction_modality:

  character of length 1, perform brain extraction modality

- template_transform_type:

  either `'Rigid'` or `'Affine'` align to template brain

- template:

  template image (not skull-stripped) or string, e.g. `'biobank'`,
  `'croppedMni152'`

- do_bias_correction:

  whether to perform bias field correction

- return_bias_field:

  return bias field as an additional output without bias correcting the
  image

- do_denoising:

  whether to remove noises using non-local means

- intensity_matching_type:

  either `'regression'` or `'histogram'`; only is performed if
  `reference_image` is not `NULL`.

- reference_image:

  `'ANTsImage'` or path to image, or `NULL`

- intensity_normalization_type:

  either re-scale the intensities to `c(0, 1)` (`'01'`), or for
  zero-mean, unit variance (`'0mean'`); if `NULL` normalization is not
  performed

- verbose:

  print progress to the screen

## Value

Dictionary with images after process. The images are registered and
re-sampled into template.

## See also

`antspynet$preprocess_brain_image`

## Examples

``` r
library(rpyANTs)
if(interactive() && ants_available("antspynet")) {
  image_path <- ants$get_ants_data('r30')
  preprocessed <- antspynet_preprocess_brain_image(
    image_path, verbose = FALSE
  )

  # Compare
  orig_img <- as_ANTsImage(image_path)
  new_img <- preprocessed$preprocessed_image
  pal <- grDevices::gray.colors(256, start = 0, end = 1)

  par(mfrow = c(1, 2), mar = c(0.1, 0.1, 0.1, 0.1),
      bg = "black", fg = "white")
  image(orig_img[], asp = 1, axes = FALSE,
        col = pal, ylim = c(1, 0))
  image(new_img[], asp = 1, axes = FALSE,
        col = pal, ylim = c(1, 0))

}

```
