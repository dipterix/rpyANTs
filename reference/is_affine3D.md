# Check if an object is a 3D 'affine' transform matrix

Check if an object is a 3D 'affine' transform matrix

## Usage

``` r
is_affine3D(x, ...)

# Default S3 method
is_affine3D(x, strict = TRUE, ...)

# S3 method for class 'ants.core.ants_transform.ANTsTransform'
is_affine3D(x, ...)
```

## Arguments

- x:

  R or Python object, accepted forms are numeric `matrix`,
  `'ANTsTransform'`, or `character` (path to transform matrix)

- ...:

  passed to other methods

- strict:

  whether the last element should be always 1

## Value

A logical value whether the object can be loaded as a 4-by-4 matrix.

## Examples

``` r
# not affine
is_affine3D(1)
#> [1] FALSE

# 3x3 matrix is not as it is treated as 2D transform
is_affine3D(matrix(rnorm(9), nrow = 3))
#> [1] FALSE

# 3x4 matrix
x <- matrix(rnorm(12), nrow = 3)
is_affine3D(x)
#> [1] TRUE

# 4x4 matrix
x <- rbind(x, c(0,0,0,1))
is_affine3D(x)
#> [1] TRUE

if(interactive() && ants_available()) {

  ants <- load_ants()
  x <- ants$new_ants_transform(dimension = 3L)
  is_affine3D(x)

  # save the parameters
  f <- tempfile(fileext = ".mat")
  ants$write_transform(x, f)
  is_affine3D(f)

}



```
