# Get `'ANTsPyNet'` module

Get `'ANTsPyNet'` module

## Usage

``` r
load_antspynet(force = FALSE, error_if_missing = TRUE)
```

## Arguments

- force:

  whether to force reloading `antspynet` module; default is false

- error_if_missing:

  whether to raise errors when the module is unable to load; default is
  true.

## Value

A 'Python' module if successfully loaded. If `error_if_missing` is set
to false and module is unable to load, return `NULL`

## See also

[`ants`](http://dipterix.org/rpyANTs/reference/ants.md)
