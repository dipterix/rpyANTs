# Resample image

See `ants$resample_image` for more details

## Usage

``` r
ants_resample_image(
  x,
  resample_params,
  use_voxels = FALSE,
  interp_type = c("linear", "nn", "guassian", "sinc", "bspline")
)
```

## Arguments

- x:

  input image

- resample_params:

  either relative number or absolute integers

- use_voxels:

  whether the `resample_params` should be treated as new dimension
  `use_voxels=TRUE`, or the new dimension should be calculated based on
  current dimension and `resample_params` combined (`use_voxels=FALSE`
  then `resample_params` will be treated as relative number); default is
  `FALSE`

- interp_type:

  interpolation type; either integer or character; see 'Usage' for
  available options

## Value

Resampled image

## Examples

``` r

if(interactive() && ants_available()) {

  sample_fpath <- as.character(ants$get_ants_data("r16"))

  if(file.exists(sample_fpath)) {

    try(silent = TRUE, {

      ants <- load_ants()
      fi <- as_ANTsImage(ants$get_ants_data("r16"))

      # linear (interp_type = 0 or "linear)
      filin <- ants_resample_image(fi, c(50, 60), TRUE, "linear")

      # nearest neighbor (interp_type = 1 or "nn)
      finn <- ants_resample_image(fi, c(50, 60), TRUE, "nn")

      par(mfrow = c(1, 3), mar = c(0, 0, 0, 0))
      pal <- gray.colors(256, start = 0)

      image(fi[], asp = 1, axes = FALSE,
            ylim = c(1,0), col = pal)
      image(filin[], asp = 1, axes = FALSE,
            ylim = c(1,0), col = pal)
      image(finn[], asp = 1, axes = FALSE,
            ylim = c(1,0), col = pal)

    })

  }

}

```
