# Slice index in 'Python' arrays

Slice index in 'Python' arrays

## Usage

``` r
py_slice(...)
```

## Arguments

- ...:

  passing to `slice` ('Python')

## Value

Index slice instance

## Examples

``` r

if(interactive() && ants_available()) {

  x <- np_array(array(seq(20), c(4, 5)))

  # equivalent to x[::2]
  x[py_slice(NULL, NULL, 2L)]

}
```
