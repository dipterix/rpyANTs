# Process 'T1' image

Process 'MRI' and align with template brains

## Usage

``` r
t1_preprocess(
  t1_path,
  templates = "mni_icbm152_nlin_asym_09a",
  work_path = ".",
  verbose = TRUE
)
```

## Arguments

- t1_path:

  path to a 'T1' image

- templates:

  template to use; default is `'mni_icbm152_nlin_asym_09a'`,

- work_path:

  working path, must be a directory

- verbose:

  whether to verbose the progress

## Value

Nothing will be returned. Please check `work_path` for results.
