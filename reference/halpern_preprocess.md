# 'ANTs' functions for 'Halpern' lab

'ANTs' functions for 'Halpern' lab

## Usage

``` r
halpern_register_ct_mri(
  fixed,
  moving,
  outprefix,
  fixed_is_ct = TRUE,
  verbose = TRUE
)

halpern_register_template_mri(
  fixed,
  moving,
  outprefix,
  mask = NULL,
  verbose = TRUE
)

halpern_apply_transform_template_mri(roi_folder, outprefix, verbose = TRUE)
```

## Arguments

- fixed:

  fixed image as template

- moving:

  moving image that is to be registered into `fixed`

- outprefix:

  output prefix, needs to be absolute path prefix

- fixed_is_ct:

  whether `fixed` is 'CT'

- verbose:

  whether to verbose the progress; default is true

- mask:

  mask file for template (skull-stripped)

- roi_folder:

  template 'ROI' or atlas folder in which the image atlases or masks
  will be transformed into subject's native brain

## Value

A list of result configurations
