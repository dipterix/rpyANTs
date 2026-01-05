# Extract brain and strip skull

Function `brain_mask` and `brain_extraction` use the `ants` base
package. Function `antspynet_brain_extraction` uses `antspynet` to
extract the brain using deep neural network. This requires additional
configuration. Print `antspynet$brain_extraction` to see the original
documentation.

## Usage

``` r
antspynet_brain_extraction(
  x,
  modality = c("t1", "t1nobrainer", "t1combined", "flair", "t2", "t2star", "bold", "fa",
    "t1t2infant", "t1infant", "t2infant"),
  verbose = FALSE
)

brain_mask(
  x,
  work_path = tempfile(pattern = "rpyant_brain_extraction_"),
  verbose = TRUE,
  auto_clean = TRUE
)

brain_extraction(
  x,
  skull_alpha = 0,
  threshold_quantile = 0.85,
  ...,
  verbose = TRUE
)
```

## Arguments

- x:

  input image or image path

- modality:

  modality type, used by `antspynet_brain_extraction` only

- verbose:

  whether to print out process to the screen

- work_path:

  working directory; default is temporary path

- auto_clean:

  whether to automatically clean the working path if the path is
  temporary; default is true

- skull_alpha:

  used by `brain_extraction`, the opacity of the skull in the final
  image; default is 0 (completely strip the skull); set to values
  between 0 to 1 to add skulls (in such case, the background noises will
  be set to 0)

- threshold_quantile:

  used only when `skull_alpha` is positive to remove the noises

- ...:

  see `work_path` and `auto_clean`

## Value

- `brain_mask`:

  Brain mask image

- `brain_extraction`:

  extracted brain

- `antspynet_brain_extraction`:

  Brain mask image
