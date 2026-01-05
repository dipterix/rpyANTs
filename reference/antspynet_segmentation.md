# Imaging segmentation using `antspynet`

Supports `Desikan-Killiany-Tourville` labeling and deep `'Atropos'`.

## Usage

``` r
antspynet_desikan_killiany_tourville_labeling(
  x,
  do_preprocessing = TRUE,
  return_probability_images = FALSE,
  do_lobar_parcellation = FALSE,
  verbose = TRUE
)

antspynet_deep_atropos(
  x,
  do_preprocessing = TRUE,
  use_spatial_priors = TRUE,
  aseg_only = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  `'NIfTI'` image or path to the image that is to be segmented

- do_preprocessing:

  whether `x` is in native space and needs the be registered to template
  brain before performing segmentation; default is true since the model
  is trained with template brain. If you want to manually process the
  image, see
  [`antspynet_preprocess_brain_image`](http://dipterix.org/rpyANTs/reference/antspynet_preprocess_brain_image.md)

- return_probability_images:

  whether to return probability images

- do_lobar_parcellation:

  whether to perform lobar 'parcellation'

- verbose:

  whether to print out the messages

- use_spatial_priors:

  whether to use `'MNI'` partial tissue priors

- aseg_only:

  whether to just return the segmented image

## Value

One or a list of `'ANTsImage'` image instances. Please print out
`antspynet$desikan_killiany_tourville_labeling` or
`antspynet$deep_atropos` to see the details.

## See also

`antspynet$desikan_killiany_tourville_labeling`,
`antspynet$deep_atropos`

## Examples

``` r

# Print Python documents
if(interactive() && ants_available("antspynet")) {
  antspynet <- load_antspynet()

  print(antspynet$deep_atropos)

  print(antspynet$desikan_killiany_tourville_labeling)
}

```
