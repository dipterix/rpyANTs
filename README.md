
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rpyANTs

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/rpyANTs)](https://CRAN.R-project.org/package=rpyANTs)
[![R-check](https://github.com/dipterix/rpyANTs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dipterix/rpyANTs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`rpyANTs` is a package that ports `ANTsPy`, a `Python` implementation of
`ANTs` into R using R-Python interpreter package `reticulate`.

> Disclaimer: This is NOT the official [`ANTsR`
> package](https://github.com/ANTsX/ANTsR).

## Installation

The installation requires one-line extra setup

``` r
# Install from CRAN
install.packages("rpyANTs")

# Install from nightly dev builder
# install.packages("rpyANTs", repos = "https://dipterix.r-universe.dev")


# set up ANTs
rpyANTs::install_ants()
```

`install_ants` creates an isolated `Python` environment managed by
[`RAVE`](https://openwetware.org/wiki/RAVE). This environment does not
conflict nor affect your existing Python installations.

## How to use

To load `ANTs`

``` r
library(rpyANTs)

# Whether ANTs is available
ants_available()

# Load ANTs into R
ants
```

In R, we use `$` to get module functions or class members. For example:

``` r
ants$add_noise_to_image
#> <ANTs Python Wrapper>
#> Help on function add_noise_to_image in module ants.utils.add_noise_to_image:
#> 
#> add_noise_to_image(image, noise_model, noise_parameters)
#>     Add noise to an image using additive Guassian, salt-and-pepper,
#>     shot, or speckle noise.
#>     
#>     ANTsR function: `addNoiseToImage`
#>     
#>     Arguments
#>     ---------
#>     image : ANTsImage
#>         scalar image.
#>     
#>     noise_model : string
#>         'additivegaussian', 'saltandpepper', 'shot', or 'speckle'.
#>     
#>     noise_parameters : tuple or array or float
#>         'additivegaussian': (mean, standardDeviation)
#>         'saltandpepper': (probability, saltValue, pepperValue)
#>         'shot': scale
#>         'speckle': standardDeviation
#>     
#>     Returns
#>     -------
#>     ANTsImage
#>     
#>     Example
#>     -------
#>     >>> import ants
#>     >>> image = ants.image_read(ants.get_ants_data('r16'))
#>     >>> noise_image = ants.add_noise_to_image(image, 'additivegaussian', (0.0, 1.0))
#>     >>> noise_image = ants.add_noise_to_image(image, 'saltandpepper', (0.1, 0.0, 100.0))
#>     >>> noise_image = ants.add_noise_to_image(image, 'shot', 1.0)
#>     >>> noise_image = ants.add_noise_to_image(image, 'speckle', 1.0)
#> 
#> *** Above documentation is for Python. 
#> *** Please use `$` instead of `.` for modules and functions in R
#> <function add_noise_to_image at 0x115d50790>
```

The following R code translates Python code into R:

``` r
# >>> img = ants.image_read(ants.get_ants_data('r16'))
img <- ants$image_read(ants$get_ants_data('r16'))

# >>> noise_image1 = ants.add_noise_to_image(img, 'additivegaussian', (0.0, 1.0))
noise_image1 <- ants$add_noise_to_image(
  img, 'additivegaussian', 
  noise_parameters = tuple(0.0, 1.0)
)

# >>> noise_image2 = ants.add_noise_to_image(img, 'saltandpepper', (0.1, 0.0, 100.0))
noise_image2 <- ants$add_noise_to_image(
  img, 'saltandpepper', 
  noise_parameters = tuple(0.1, 0.0, 100.0)
)

# >>> noise_image3 = ants.add_noise_to_image(img, 'shot', 1.0)
noise_image3 <- ants$add_noise_to_image(
  img, 'shot', 
  noise_parameters = 1.0
)

# >>> noise_image4 = ants.add_noise_to_image(img, 'speckle', 1.0)
noise_image4 <- ants$add_noise_to_image(
  img, 'speckle', 
  noise_parameters = 1.0
)
```

To load imaging data into R

``` r
orig_array <- py_to_r(img$numpy())
noise_array1 <- py_to_r(noise_image1$numpy())
noise_array2 <- py_to_r(noise_image2$numpy())
noise_array3 <- py_to_r(noise_image3$numpy())
noise_array4 <- py_to_r(noise_image4$numpy())

# plot via R
layout(matrix(c(1,1,2,3,1,1,4,5), nrow = 2, byrow = TRUE))
par(mar = c(0.1, 0.1, 0.1, 0.1), bg = "black", fg = "white")
pal <- grDevices::gray.colors(256, start = 0, end = 1)

image(orig_array, asp = 1, axes = FALSE, 
      col = pal, zlim = c(0, 255), ylim = c(1, 0))
image(noise_array1, asp = 1, axes = FALSE, 
      col = pal, zlim = c(0, 255), ylim = c(1, 0))
image(noise_array2, asp = 1, axes = FALSE, 
      col = pal, zlim = c(0, 255), ylim = c(1, 0))
image(noise_array3, asp = 1, axes = FALSE, 
      col = pal, zlim = c(0, 255), ylim = c(1, 0))
image(noise_array4, asp = 1, axes = FALSE, 
      col = pal, zlim = c(0, 255), ylim = c(1, 0))
```

<img src="man/figures/README-load-image-into-r-1.png" width="80%" />

## Upgrade `ANTs`

To upgrade `ANTs`, first update `rpyANTs`, then upgrade `ANTsPyx`

``` r
install.packages("rpyANTs")
rpymat::add_packages(packages = "antspyx", pip = TRUE)
```

## Advanced use case

#### Run/Debug `Python` scripts

`rpyANTs` ports functions that allows to run `Python` scripts. For
example:

``` r
library(rpyANTs)

script_path <- tempfile(fileext = ".py")
writeLines(con = script_path, text = r"(

# This is Python script
import ants
print(ants.__version__)

)")

run_script(script_path)
```

You can also run `Python` interactive in R (yes, you are correct).
Simply run

``` r
rpyANTs::repl_python()
```

The console prefix will change from `>` to `>>>`, meaning you are in
`Python` mode:

    > rpyANTs::repl_python()
    Python 3.8.16 (/Users/dipterix/Library/r-rpymat/miniconda/envs/rpymat-conda-env/bin/python3.8)
    Reticulate 1.26 REPL -- A Python interpreter in R.
    Enter 'exit' or 'quit' to exit the REPL and return to R.
    >>> 

Try some Python code!

    >>> import ants
    >>> help(ants.registration)

To exit Python mode, type `exit` (no parenthesis) and hit enter key

    >>> exit
    > 

#### Data conversions

Native R variables can be easily converted to `Python` and back via
`r_to_py` and `py_to_r`.

For example

``` r
# R to Python
r_to_py(1)
#> 1.0
r_to_py(1L)
#> 1

# Python to R
py_obj <- py_list(1:3)
class(py_obj)  # <- this is a python object
#> [1] "python.builtin.list"   "python.builtin.object"

py_to_r(py_obj)
#> [1] 1 2 3
```

You can also use variables created in R from Python or vice versa:

In the following example, an R object `object_r` is created. In Python,
it can be accessed (read-only) via `r.object_r`

    > object_r <- c(1,2,3)
    > repl_python()
    Python 3.8.16 (/Users/dipterix/Library/r-rpymat/miniconda/envs/rpymat-conda-env/bin/python3.8)
    Reticulate 1.26 REPL -- A Python interpreter in R.
    Enter 'exit' or 'quit' to exit the REPL and return to R.
    >>> r.object_r
    [1.0, 2.0, 3.0]

Similarly, a Python object `object_py` is created, and it can be read
from `py$object_py`:

    >>> import numpy as np
    >>> object_py = np.array([2,3,4])
    >>> exit
    > py$object_py
    [1] 2 3 4

## Known issues

#### Variable types

R is not a type-rigid language. Some functions in `ANTsPy` require
specific variable types that are often vague in R. For example the
`dimension` argument in function `ants$create_ants_transform` needs to
be an integer, but R’s default numerical values are `double`. In this
case, variable formats need to be explicitly given.

Here are several examples

1.  Explicit integers

``` r
# ants$create_ants_transform(dimension = 3)     # <- error
ants$create_ants_transform(dimension = 3L)      # < XXXL is an explicit integer
```

2.  `Tuple`, `list`, and `dictionary`

A Python `tuple` is a vector that cannot alter lengths.

``` r
# Wrong as `aff_iterations` needs to be a tuple
# ants$registration(fixed, moving, ..., aff_iterations = c(6L, 4L, 2L, 1L))

ants$registration(fixed, moving, ..., aff_iterations = tuple(6L, 4L, 2L, 1L))
```

Similar conversions can be done via `py_list`, `py_dict`.

#### Operators

In Python, operators on `ANTsImage`, such as `img > 5` are defined. Such
operators is being supported in R as `S3` generic functions. Don’t worry
if you don’t know what is `S3` generic, see the following examples:

``` r
library(rpyANTs)
image <- ants$image_read(ants$get_ants_data('mni'))
print(image)
dim(image)
range(image)

y1 <- (image > 10) * 8000

y2 <- image
y2[y2 < 10] <- 4000

y3 <- log(image + 1000)
y3 <- (y3 - min(y3)) / (max(y3) - min(y3)) * 8000

ants_plot_grid(
  list(image, y1, y2, y3),
  slices = 100, shape = c(2,2),
  vmin = 0, vmax = 8000
)
```

<img src="man/figures/README-s3-generic-showcase-1.png" width="80%" />

Notice the operator generics are still under implementation. Some
classes/objects are still not supported. In this case, you might want to
use the following workaround methods. You are more than welcome to post
a wish-list or issue ticket to the [`Github`
repository](https://github.com/dipterix/rpyANTs/issues)

Alternative version 1: call operators directly

``` r
library(rpyANTs)
image <- ants$image_read(ants$get_ants_data('r16'))

# The followings are the same
# threshold <- image > 10
threshold <- image$`__gt__`(10)
ants$plot(threshold)
```

Work-around version 2: If you don’t know how Python operators work, use
Python directly

``` r
library(rpyANTs)
image <- ants$image_read(ants$get_ants_data('r16'))

# Create an R variable from Python!
py_run_string("r.threshold = r.image > 10", local = TRUE, convert = FALSE)
ants$plot(threshold)
```

## Citation

This is a general citation:

> Avants, B.B., Tustison, N. and Song, G., 2009. Advanced normalization
> tools (ANTS). The Insight Journal, 2(365), pp.1-35.

Please check their official website to cite specific methods.

## License

This package `rpyANTs` is released under Apache-2.0 license (Copyright:
Zhengjia Wang). The underlying `ANTsPy` is released under Apache-2.0
license (Copyright: ANTs contributors).
