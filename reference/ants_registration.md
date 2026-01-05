# Register two images using 'ANTs'

Register two images using 'ANTs'

## Usage

``` r
ants_registration(
  fixed,
  moving,
  type_of_transform = "SyN",
  initial_transform = NULL,
  outprefix = tempfile(),
  mask = NULL,
  grad_step = 0.2,
  flow_sigma = 3,
  total_sigma = 0,
  aff_metric = c("mattes", "GC", "meansquares"),
  aff_sampling = 32,
  aff_random_sampling_rate = 0.2,
  syn_metric = c("mattes", "CC", "meansquares", "demons"),
  syn_sampling = 32,
  reg_iterations = c(40, 20, 0),
  aff_iterations = c(2100, 1200, 1200, 10),
  aff_shrink_factors = c(6, 4, 2, 1),
  aff_smoothing_sigmas = c(3, 2, 1, 0),
  write_composite_transform = FALSE,
  verbose = FALSE,
  smoothing_in_mm = FALSE,
  ...
)
```

## Arguments

- fixed:

  fixed image to which we register the moving image, can be character
  path to 'NIfTI' image, or `'ANTsImage'` instance, `'oro.nifti'`
  object, `'niftiImage'` from package `'RNifti'`, or `'threeBrain.nii'`
  from package `'threeBrain'`; see also
  [`as_ANTsImage`](http://dipterix.org/rpyANTs/reference/as_ANTsImage.md)

- moving:

  moving image to be mapped to fixed space; see also
  [`as_ANTsImage`](http://dipterix.org/rpyANTs/reference/as_ANTsImage.md)

- type_of_transform:

  a linear or non-linear registration type; print `ants$registration` to
  see details

- initial_transform:

  optional list of strings; transforms to apply prior to registration

- outprefix:

  output file to save results

- mask:

  image mask; see also
  [`as_ANTsImage`](http://dipterix.org/rpyANTs/reference/as_ANTsImage.md)

- grad_step, flow_sigma, total_sigma:

  optimization parameters

- aff_metric:

  the metric for the 'affine' transformation, choices are `'GC'`,
  `'mattes'`, `'meansquares'`

- aff_sampling, aff_random_sampling_rate, aff_iterations,
  aff_shrink_factors, aff_smoothing_sigmas:

  controls 'affine' transform

- syn_metric:

  the metric for the `'SyN'` transformation, choices are `'GC'`,
  `'mattes'`, `'meansquares'`, `'demons'`

- syn_sampling, reg_iterations:

  controls the `'SyN'` transform

- write_composite_transform:

  whether the composite transform (and its inverse, if it exists) should
  be written to an 'HDF5' composite file; default is false

- verbose:

  verbose the progress

- smoothing_in_mm:

  logical, currently only impacts low dimensional registration

- ...:

  others passed to `ants$registration`

## Value

A 'Python' dictionary of aligned images and transform files.

## Details

Function family `ants_registration*` align images (specified by
`moving`) to `fixed`. Here are descriptions of the variations:

- `ants_registration`:

  Simple wrapper function for 'Python' implementation
  `ants.registration`, providing various of registration options

- `ants_registration_halpern1`:

  Rigid-body registration designed for 'Casey-Halpern' lab, mainly used
  for aligning 'MRI' to 'CT' (or the other way around)

## Examples

``` r
if(interactive() && ants_available()) {

  ants <- load_ants()

  # check the python documentation here for detailed explanation
  print(ants$registration)

  # example to register
  fi <- ants$image_read(ants$get_ants_data('r16'))
  mo <- ants$image_read(ants$get_ants_data('r64'))

  # resample to speed up this example
  fi <- ants$resample_image(fi, list(60L,60L), TRUE, 0L)
  mo <- ants$resample_image(mo, list(60L,60L), TRUE, 0L)

  # SDR transform
  transform <- ants_registration(
    fixed=fi, moving=mo, type_of_transform = 'SyN' )

  ants$plot(fi, overlay = transform$warpedmovout, overlay_alpha = 0.3)


}

```
