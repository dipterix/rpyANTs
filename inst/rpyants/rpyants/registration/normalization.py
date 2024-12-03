#########################################################
# Custom Normalization scripts
#########################################################

import ants
from ..utils.paths import normalize_path
from ..utils.internals import get_lib_fn, ants_process_arguments

# from rpyants.registration.normalization import normalize_to_template_syn
# fix_path="/Users/dipterix/rave_data/raw_dir/PAV044/rave-imaging/normalization/anat/MNI152NLin2009bAsym_orig.nii.gz"
# moving_path="/Users/dipterix/rave_data/raw_dir/Precision001/rave-imaging/inputs/MRI/MRI_RAW.nii"
# outprefix="/Users/dipterix/rave_data/raw_dir/Precision001/rave-imaging/mornalize_mni152b_"
# normalize_to_template_syn(fix_path, moving_path, outprefix, verbose=True)

def normalize_to_template_syn(
    fix_path : str, mov_paths : list | tuple, outprefix : str, 
    weights : float | int | list | tuple = 1, verbose : bool = True,
    rigid_parameters = {
      "convergence" : "[1000x500x250x0,1e-7,10]",
      "shrink-factors" : "12x8x4x1",
      "smoothing-sigmas" : "5x4x3x1vox",
    },
    affine_parameters = {
      "convergence" : "[1000x500x250x0,1e-7,10]",
      "shrink-factors" : "12x8x4x1",
      "smoothing-sigmas" : "5x4x3x1vox",
    },
    syn_parameters = {
      "convergence" : "[1000x500x500x0,1e-7,7]",
      "shrink-factors" : "8x4x4x1",
      "smoothing-sigmas" : "4x3x1x1vox",
    },
    refine_parameters = {
      "convergence" : "[200x50x10x0,1e-6,7]",
      "shrink-factors" : "4x4x2x1",
      "smoothing-sigmas" : "2x2x1x1vox",
    }):
    '''
    Normalize images to template via volumetric mapping in 4 satges
    
    @param fix_path: Absolute path to the fixing (template) image
    @type fix_path: str
    
    @param fix_path: A list or tuple of absolute paths to the moving (native) images; all images should be registered to the same space
    @type fix_path: list | tuple
    
    @param outprefix: Output path prefix
    @type outprefix: str
    
    @param weights: A single number or a list of numbers indicating the weights for each moving image
    @type weights: float | int | list | tuple
    
    @param rigid_parameters: dictionary of stage-1 rigid parameters
    @type rigid_parameters: dict; default is 
    {
      "convergence" : "[1000x500x250x0,1e-7,10]",
      "shrink-factors" : "12x8x4x1",
      "smoothing-sigmas" : "5x4x3x1vox",
    }
    
    @param affine_parameters: dictionary of stage-2 affine parameters
    @type affine_parameters: dict; default is
    {
      "convergence" : "[1000x500x250x0,1e-7,10]",
      "shrink-factors" : "12x8x4x1",
      "smoothing-sigmas" : "5x4x3x1vox",
    }
    
    @param syn_parameters: dictionary of stage-3 SyN parameters
    @type syn_parameters: dict; default is
    {
      "convergence" : "[1000x500x500x0,1e-7,7]",
      "shrink-factors" : "8x4x4x1",
      "smoothing-sigmas" : "4x3x1x1vox",
    }
    
    @param refine_parameters: dictionary of stage-4 refine (also SyN) parameters
    @type refine_parameters: dict; default is
    {
      "convergence" : "[200x50x10x0,1e-6,7]",
      "shrink-factors" : "4x4x2x1",
      "smoothing-sigmas" : "2x2x1x1vox",
    }
    
    @param verbose: Whether to verbose outputs
    @type verbose: bool
    '''
    # use forward slash to avoid escaping issues on windows
    # Do not use memory pointers!!!
    fix_path = normalize_path(fix_path, sep="/")
    if not isinstance(mov_paths, list):
        raise TypeError("normalize_to_template_syn: `mov_paths` must be a list")
    mov_paths = [normalize_path(moving_path, sep="/") for moving_path in mov_paths]
    outprefix = normalize_path(outprefix, sep="/")
    verbose = "1" if verbose else "0"
    if not isinstance(weights, list):
        weights = [weights for x in mov_paths]
    if len(weights) != len(mov_paths):
        raise ValueError("normalize_to_template_syn: weights must be a single number or a list of the same length as `mov_paths`")
    libfn = get_lib_fn("antsRegistration")
    # Build arguments, see https://github.com/ANTsX/ANTsPy/blob/78dad33ca4e12ae605e80d8d990ec73a52abe273/ants/registration/interface.py
    metrics = []
    for moving_path, weight in zip(mov_paths, weights):
        metrics.append("--metric")
        metrics.append("MI[%s,%s,%.6f,32,Random,0.25]" % (fix_path, moving_path, weight))
    args = [
        "--dimensionality", "3",
        "--float", "1",
        "-o", "[%s,%s,%s]" % (outprefix, outprefix + "Warped.nii.gz", outprefix + "InverseWarped.nii.gz"),
        "--interpolation", "Linear",
        "--use-histogram-matching", "1",
        "-v", verbose,
        # Stage 1: Rigid Transformation
        # --transform Rigid[0.1] \
        # --metric MI[fixedImage, movingImage, 1, 32, Random, 0.25] \
        # --convergence [1000x500x250x0,1e-7,10] \
        # --shrink-factors 12x8x4x1 \
        # --smoothing-sigmas 5x4x3x1vox \
        "--transform", "Rigid[2.0]",
    ] + metrics + [
        "--convergence", rigid_parameters.get("convergence", "[1000x500x250x0,1e-7,10]"),
        "--shrink-factors", rigid_parameters.get("shrink-factors", "12x8x4x1"),
        "--smoothing-sigmas", rigid_parameters.get("smoothing-sigmas", "5x4x3x1vox"),
        # Stage 2: Affine Transformation
        # --transform Affine[0.1] \
        # --metric MI[/path/to/MNI152b.nii.gz,/path/to/subject_T1w.nii.gz,1,32,Random,0.25] \
        # --convergence [1000x500x250x0,1e-7,10] \
        # --shrink-factors 12x8x4x1 \
        # --smoothing-sigmas 5x4x3x1vox \
        "--transform", "Affine[0.1]",
    ] + metrics + [
        "--convergence", affine_parameters.get("convergence", "[1000x500x250x0,1e-7,10]"),
        "--shrink-factors", affine_parameters.get("shrink-factors", "12x8x4x1"),
        "--smoothing-sigmas", affine_parameters.get("smoothing-sigmas", "5x4x3x1vox"),
        # Stage 3: SyN Transformation
        # --transform SyN[0.3,4,3] \
        # --metric MI[/path/to/MNI152b.nii.gz,/path/to/subject_T1w.nii.gz,1,32,Random,0.25] \
        # --convergence [1000x500x500x0,1e-7,7] \
        # --shrink-factors 8x4x4x1 \
        # --smoothing-sigmas 4x3x1x1vox
        "--transform", "SyN[0.3,4.0,3.0]",
    ] + metrics + [
        "--convergence", syn_parameters.get("convergence", "[1000x500x500x0,1e-7,7]"),
        "--shrink-factors", syn_parameters.get("shrink-factors", "8x4x4x1"),
        "--smoothing-sigmas", syn_parameters.get("smoothing-sigmas", "4x3x1x1vox"),
        # Stage 4: subcortical refine stage
        "--transform", "SyN[0.3,4,3]",
    ] + metrics + [
        "--convergence", refine_parameters.get("convergence", "[200x50x10x0,1e-6,7]"),
        "--shrink-factors", refine_parameters.get("shrink-factors", "4x4x2x1"),
        "--smoothing-sigmas", refine_parameters.get("smoothing-sigmas", "2x2x1x1vox"),
    ]
    # Clean args
    processed_args = ants_process_arguments(args)
    reg_exit = libfn(processed_args)
    if verbose == "1":
        print("---- Call arguments --------------------------------")
        print("antsRegistration.sh %s" % " ".join(processed_args))
    if not reg_exit == 0:
        print(f"Registration failed with error code {reg_exit}")
    return {
        "outprefix"    : outprefix,
        "warpedmovout" : outprefix + "Warped.nii.gz",
        "warpedfixout" : outprefix + "InverseWarped.nii.gz",
        "fwdtransforms": [
            outprefix + "1Warp.nii.gz",
            outprefix + "0GenericAffine.mat", 
        ],
        "invtransforms": [
            outprefix + "0GenericAffine.mat", 
            outprefix + "1InverseWarp.nii.gz"
        ]
    }

# Some obsolete code from yae.py

# # step 1: align moving image to fixing image with affine registration
# fix_img = ants.image_read(fix_path)
# mov_img = ants.image_read(mov_path)
# affine_reg = ants.registration(
#     fixed = fix_img, moving = mov_img, 
#     mask = template_mask, moving_mask = native_mask,
#     outprefix = file_path(ants_outputdir, "affine_"), 
#     type_of_transform = "antsRegistrationSyN[a]", 
#     write_composite_transform = False, 
#     aff_metric = "mattes",  
#     aff_random_sampling_rate=0.25,
#     aff_iterations=(1000, 500, 250, 0),
#     aff_shrink_factors=(12, 8, 4, 1),
#     aff_smoothing_sigmas=(5, 4, 3, 1),
#     verbose = verbose
# )
# # self.expected_image_path(native_type, "normalization/anat", name = "CT_resampled", space = "T1RAS")
# fix_img.to_file(ensure_basename(file_path(resultdir, "anat", f"{template_name}_orig.nii.gz")))
# # fix_img is the template
# # mov_img is the native image
# # affine_reg['warpedmovout'] is the registered T1w image in template
# affine_reg['warpedmovout'].to_file(ensure_basename(
#     self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "affine")
# ))
# affine_reg['warpedfixout'].to_file(ensure_basename(
#     self.expected_image_path(native_type, "normalization/anat", name = template_name, space = native_type, transform = "affine")
# ))
# # initial_transform = ants.read_transform(affine_reg['fwdtransforms'][0])
# multivariate_extras=[]
# for modality, path in other_images.items():
#     # ("nameOfMetric2", img, img, weight, metricParam )
#     # other_image = ants.image_read(path)
#     extra_criteria = (
#         "MI", fix_img, 
#         # initial_transform.apply_to_image(other_image),
#         ants.image_read(path),
#         self.normalization_weights.get(modality, 0.6),
#         "32,Random,0.25"
#     )
#     multivariate_extras.append(extra_criteria)
# if len(multivariate_extras) == 0:
#     multivariate_extras = None
# # Step 2: SyN normalization
# # aff_* parameters will not be used since this is SyN-only registration
# # with quick iterations
# # deformable_transform_type = "antsRegistrationSyNQuick[so]"
# deformable_transform_type = "SyNOnly"
# deformable_reg = ants.registration(
#     fixed = fix_img, moving = mov_img,
#     mask = template_mask, moving_mask = native_mask,
#     outprefix = file_path(ants_outputdir, "deformable_"), 
#     type_of_transform = deformable_transform_type,
#     flow_sigma = 4,
#     total_sigma = 3,
#     initial_transform = affine_reg['fwdtransforms'][0],
#     write_composite_transform = False,
#     aff_metric = "mattes", 
#     aff_random_sampling_rate=0.25,
#     aff_iterations=(1000, 500, 250, 0),
#     aff_shrink_factors=(12, 8, 4, 1),
#     aff_smoothing_sigmas=(5, 4, 3, 1),
#     syn_metric="mattes", 
#     reg_iterations=(1000, 600, 500, 10, 0),
#     multivariate_extras = multivariate_extras,
#     verbose = verbose
# )
# # save forward (mov -> fix) transforms
# # simple experiment to check the composition of the transforms: the following two are the same
# # np.matmul(ants_AffineTransform_to_m44(affine_1), ants_AffineTransform_to_m44(affine_0))
# # ants_AffineTransform_to_m44(ants.compose_ants_transforms([affine_0, affine_1]))
# file_copy( affine_reg['fwdtransforms'][0], ensure_basename(
#     file_path( self._work_path, "normalization/transformations", transform_to_template_prefix + "ants0.mat" )
# ))
# file_copy( deformable_reg['fwdtransforms'][0], ensure_basename(
#     file_path( self._work_path, "normalization/transformations", transform_to_template_prefix + "ants1.nii.gz" )
# ))
# # file_copy( deformable_reg['fwdtransforms'][1], ensure_basename(
# #     file_path( self._work_path, "normalization/transformations", transform_to_template_prefix + "ants1.mat" )
# # ))
# # file_copy( deformable_reg['fwdtransforms'][0], ensure_basename(
# #     file_path( self._work_path, "normalization/transformations", transform_to_template_prefix + "ants2.nii.gz" )
# # ))
# # calculate and save the inverse transforms (fix -> mov)
# file_copy( deformable_reg['invtransforms'][1], ensure_basename(
#     file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants0.nii.gz" )
# ))
# ants.write_transform(
#     ants.read_transform(deformable_reg['fwdtransforms'][1]).invert(),
#     file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants1.mat" )
# )
# # ants.write_transform(
# #     ants.read_transform(affine_reg['fwdtransforms'][0]).invert(),
# #     file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants2.mat" )
# # )
# # help(ants.apply_transforms)
# # The commented code reproduces the same result as `deformable_reg['warpedmovout']`
# # warped_T1_img = ants.apply_transforms(
# #     fixed = fix_img, moving = mov_img,
# #     transformlist=[deformable_reg['fwdtransforms'][0], deformable_reg['fwdtransforms'][1], affine_reg['fwdtransforms'][0]],
# #     whichtoinvert = [False, False, False]
# # )
# # warped_T1_img.to_file(ensure_basename(
# #     self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "deformable")
# # ))
# deformable_reg['warpedmovout'].to_file(ensure_basename(
#     self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "deformable")
# ))
# # fix_img.plot(overlay = deformable_reg['warpedmovout'].iMath_canny(0.1, -1, 1), overlay_alpha = 0.5)
# warped_template_img = ants.apply_transforms(
#     fixed = mov_img, moving = fix_img,
#     transformlist=[
#         # file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants2.mat" ),
#         file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants1.mat" ),
#         file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants0.nii.gz" )
#     ],
#     whichtoinvert = [
#         # False,
#         False, 
#         False
#     ]
# )
# warped_template_img.to_file(ensure_basename(
#     self.expected_image_path(native_type, "normalization/anat", name = template_name, space = native_type, transform = "deformable")
# ))
# mappings = self.get_template_mapping(template_name = template_name, native_type = native_type, relative = True)
# if mappings is not None:
#     log_file = ensure_basename( self.format_path(folder="normalization/log", name="mappings", ext="json", 
#                                                 template = template_name, native = native_type) )
#     with open(log_file, "w") as f:
#         json.dump(mappings, f)
# return mappings
