# Get 'ANTsPy' module

Get 'ANTsPy' module

## Usage

``` r
ants

load_ants(force = FALSE, error_if_missing = TRUE)
```

## Arguments

- force:

  whether to force reloading `ants` module; default is false

- error_if_missing:

  whether to raise errors when the module is unable to load; default is
  true.

## Value

A 'Python' module if successfully loaded. If `error_if_missing` is set
to false and module is unable to load, return `NULL`

## See also

[`antspynet`](http://dipterix.org/rpyANTs/reference/antspynet.md)
