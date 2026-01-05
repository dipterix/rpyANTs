# List in 'Python'

List in 'Python'

## Usage

``` r
py_list(..., convert = FALSE)
```

## Arguments

- ...:

  passing to `list` ('Python')

- convert:

  whether to convert the results back into R; default is no

## Value

List instance, or an R vector if converted

## Examples

``` r

if(interactive() && ants_available()) {

  py_list(list(1,2,3))
  py_list(c(1,2,3))

  py_list(array(1:9, c(3,3)))
  py_list(list(list(1:3), letters[1:3]))

}
```
