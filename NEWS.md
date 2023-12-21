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
