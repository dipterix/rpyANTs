# rpyANTs 0.0.5

* Added backward compatibility support for `Py39`
* Fixed the missing `pow` function from `numpy` 
* Added `atropos` segments for templates
* New normalization: added `atropos` segments for normalization to templates
* `ANTsPyNet` model locations are stored in package data directory to comply with `CRAN` policy
* Always run `abp` during normalization
* Added brain extraction function using base `ANTs`
* Fixed issue when subject code contains underscore, causing non-comply `BIDS` format, causing parsing errors
* Added `fsaverage` template for normalization

# rpyANTs 0.0.4

* Added `get_xform` and `set_xform` to get and set `qform` in `ANTsImage`
* Removed dependence `nibabel`, which might conflict with `ANTsPy` on certain platforms
* Using the full path to `ANTsImage` class instead of the shortcut one from `Python` to avoid issues where `ANTsImage` is not exported (`ANTsPyx==0.5.3`)
* Using multiple imaging files and faster normalization parameters when normalizing the native images to template
* Fixed `get_lib_fn` not available under `ANTsPyx>=0.5.3`
* Changed default `Python` version
* Fixed bugs when mapping native `ROI` to template
* Supported `YAEL` processing pipeline and use `ANTsImage` as inputs
* Added mapping images from native to template
* Fixed a typo in folder naming
* Method `get_native_mapping` gets correct absolute paths
* Allowed `YAEL` class to add default images through constructor
* Added intensity correction
* Exported `ensure_template` and added extra information
* Set time-out to ensure downloading the correct templates
* Added `YAEL` class to support pipelines mentioned in the paper for image registration and normalization

# rpyANTs 0.0.3

* Added registration code that reproduces `Halpern` lab results
* Added `t1_preprocess` to align `MRI` to template using non-linear registration
* No more injection to `Python` on `Windows` as the memory pointer crash has been fixed by upstream developers
* No active binding to `antspynet` as this `Python` package might fail to load on some systems

# rpyANTs 0.0.2

* Added `antspynet` to the installation
* Added `ants_apply_transforms` and `ants_apply_transforms_to_points` to apply registration transforms
* Added `ants_resample_image` to resample images
* Added `antspynet_brain_extraction` to extract brain and strip skulls
* Added `antspynet_deep_atropos` and  `antspynet_desikan_killiany_tourville_labeling` to segment brain images

# rpyANTs 0.0.1

* Added a `NEWS.md` file to track changes to the package.
