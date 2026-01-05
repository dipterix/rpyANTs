# Get 'Python' built-in object

Get 'Python' built-in object

## Usage

``` r
py_builtin(name, convert = TRUE)
```

## Arguments

- name:

  object name

- convert:

  see
  [`import_builtins`](https://rstudio.github.io/reticulate/reference/import.html)

## Value

A python built-in object specified by `name`

## Examples

``` r
if(interactive() && ants_available()) {


# ------ Basic case: use python `int` as an R function ---------
py_int <- py_builtin("int")

# a is an R object now
a <- py_int(9)
print(a)
class(a)

# ------ Use python `int` as a Python function -----------------
py_int2 <- py_builtin("int", convert = FALSE)

# b in a python object
b <- py_int2(9)

# There is no '[1] ' when printing
print(b)
class(b)

# convert to R object
py_to_r(b)



}
```
