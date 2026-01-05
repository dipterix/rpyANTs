# Load data as `'ANTsImage'` class

Load data as `'ANTsImage'` class

## Usage

``` r
as_ANTsImage(x, strict = FALSE)
```

## Arguments

- x:

  data to be converted; this can be an `'ANTsImage'` instance,
  character, `'oro.nifti'` object, `'niftiImage'` from package
  `'RNifti'`, or `'threeBrain.nii'` from package `'threeBrain'`

- strict:

  whether `x` should not be `NULL`

## Value

An `'ANTsImage'` instance; use `ants$ANTsImage` to see the 'Python'
documentation

## Examples

``` r
if(interactive() && ants_available()) {

  ants <- load_ants()

  # Python string
  x1 <- ants$get_ants_data('r16')
  as_ANTsImage( x1 )

  # R character
  nii_path <- system.file(package = "RNifti",
                          "extdata", "example.nii.gz")
  as_ANTsImage( nii_path )

  # niftiImage object
  x2 <- RNifti::readNifti(nii_path)
  as_ANTsImage( x2 )

}
```
