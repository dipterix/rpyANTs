# Convert to `'ANTsTransform'`

Convert to `'ANTsTransform'`

## Usage

``` r
as_ANTsTransform(x, ...)

# Default S3 method
as_ANTsTransform(x, dimension = 3, ...)

# S3 method for class 'ants.core.ants_transform.ANTsTransform'
as_ANTsTransform(x, ...)

# S3 method for class 'ants.core.ants_image.ANTsImage'
as_ANTsTransform(x, ...)

# S3 method for class 'numpy.ndarray'
as_ANTsTransform(x, ...)

# S3 method for class 'character'
as_ANTsTransform(x, ...)
```

## Arguments

- x:

  'affine' matrix or `'numpy'` array, character path to the matrix,
  `'ANTsTransform'`, `'ANTsImage'` as displacement field.

- ...:

  passed to other methods

- dimension:

  expected transform space dimension; default is 3

## Value

An `'ANTsTransform'` object

## Examples

``` r
if(interactive() && ants_available()) {

  mat <- matrix(c(
    0, -1, 0, 128,
    1, 0, 0, -128,
    0, 0, -1, 128,
    0, 0,  0,   1
  ), ncol = 4, byrow = TRUE)

  trans <- as_ANTsTransform(mat)
  trans

  # apply transform
  trans$apply_to_point(c(120, 400, 1))

  # same results
  mat %*% c(120, 400, 1, 1)

  trans[] == mat

}
```
