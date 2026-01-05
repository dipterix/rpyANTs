# Motion correction

Print `ants$motion_correction` to see the original document

## Usage

``` r
ants_motion_correction(
  x,
  fixed = NULL,
  type_of_transform = "BOLDRigid",
  mask = NULL,
  fdOffset = 50,
  outprefix = "",
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  input image, usually 'fMRI' series

- fixed:

  fixed image to register all timepoints to

- type_of_transform:

  see
  [`ants_registration`](http://dipterix.org/rpyANTs/reference/ants_registration.md)

- mask:

  mask for image

- fdOffset:

  offset value to use in frame-wise displacement calculation

- outprefix:

  save path

- verbose:

  whether to verbose the messages

- ...:

  passed to registration methods

## Value

Motion-corrected image

## Examples

``` r
if(interactive() && ants_available()) {
  fi <- as_ANTsImage(ants$get_ants_data('ch2'))
  mytx <- ants_motion_correction( fi )

  par(mfrow = c(1, 2), mar = c(1,1,1,1))
  image(fi[,,91], asp = 1, axes = FALSE)
  image(mytx$motion_corrected[,,91], asp = 1, axes = FALSE)
}
```
