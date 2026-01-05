# Extract transform in right-anterior-superior coordinate system

Extract transform from 'ANTs' and moving image, save the transform
matrices.

## Usage

``` r
extract_coregistration(transform_path, moving_img, outprefix = NULL)
```

## Arguments

- transform_path:

  path to the transform (`'ANTs_0GenericAffine.mat'` matrix) file

- moving_img:

  path to the moving image, for example, if `'CT'` is aligned to
  `'T1w'`, then `'CT'` is the moving image

- outprefix:

  output prefix, where the matrices are saved to; default is at the same
  folder as `transform_path`

## Value

A list of two matrices: `new_sform` is the new index to image transform
(see `'sform'` in 'NIfTI' header definition). This is a four-by-four
matrix transforming from the voxel index (moving image) to the fixed
image `'RAS'` coordinates. `ras_transform` is a four-by-four matrix from
the moving image to the fixed image in `'RAS'` coordinates.
