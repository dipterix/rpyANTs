# Truncate and correct 'MRI' intensity

Uses `ants.abp_n4` to truncate and correct intensity

## Usage

``` r
correct_intensity(image, mask = NULL, intensity_truncation = c(0.025, 0.975))
```

## Arguments

- image:

  'MRI' image to be corrected, will be passed to
  [`as_ANTsImage`](http://dipterix.org/rpyANTs/reference/as_ANTsImage.md)

- mask:

  binary mask image

- intensity_truncation:

  numerical length of two, quantile probabilities to truncate.

## Value

An `'ANTsImage'` instance

## Examples

``` r
if(interactive() && ants_available()) {
  ants <- load_ants()
  scale <- (0.1 + outer(
    seq(0, 1, length.out = 256)^6,
    seq(0, 1, length.out = 256)^2,
    FUN = "+"
  )) / 6
  img = ants$image_read(ants$get_ants_data('r16')) * scale

  corrected <- correct_intensity(img)

  pal <- gray.colors(255, start = 0)
  par(mfrow = c(1, 2), mar = c(0.1, 0.1, 2.1, 0.1),
      bg = "black", fg = "white")
  image(img[], asp = 1, axes = FALSE,
        col = pal, ylim = c(1, 0),
        main = "Original", col.main = "white")
  image(corrected[], asp = 1, axes = FALSE,
        col = pal, ylim = c(1, 0),
        main = "Corrected", col.main = "white")

}
```
