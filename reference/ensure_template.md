# Ensure the template directory is downloaded

Ensure the template directory is downloaded

## Usage

``` r
ensure_template(name = BUILTIN_TEMPLATES)
```

## Arguments

- name:

  name of the template, commonly known as `'MNI152'` templates; choices
  are `"mni_icbm152_nlin_asym_09a"`, `"mni_icbm152_nlin_asym_09b"`, and
  `"mni_icbm152_nlin_asym_09c"`.

## Value

The downloaded template path

## Examples

``` r
# Do not run for testing as this will download the template
if(FALSE) {

# Default is `mni_icbm152_nlin_asym_09a`
ensure_template()
ensure_template("mni_icbm152_nlin_asym_09a")

# Using MNI152b
ensure_template("mni_icbm152_nlin_asym_09b")

}

```
