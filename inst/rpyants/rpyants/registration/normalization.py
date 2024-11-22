#########################################################
# Custom Normalization scripts
#########################################################

import ants
from ..utils.paths import normalize_path
from ..utils.internals import get_lib_fn, ants_process_arguments

# from rpyants.registration.normalization import normalize_to_template_syn
# fixed_path="/Users/dipterix/rave_data/raw_dir/PAV044/rave-imaging/normalization/anat/MNI152NLin2009bAsym_orig.nii.gz"
# moving_path="/Users/dipterix/rave_data/raw_dir/Precision001/rave-imaging/inputs/MRI/MRI_RAW.nii"
# outprefix="/Users/dipterix/rave_data/raw_dir/Precision001/rave-imaging/mornalize_mni152b_"
# normalize_to_template_syn(fixed_path, moving_path, outprefix, verbose=True)

def normalize_to_template_syn(fixed_path, moving_path, outprefix, verbose=True):
  '''
  Ridig-body coregistration for CT and MRI for Casey Halpern's lab.
  '''
  # use forward slash to avoid escaping issues on windows
  # Do not use memory pointers!!!
  fixed_path = normalize_path(fixed_path, sep="/")
  moving_path = normalize_path(moving_path, sep="/")
  outprefix = normalize_path(outprefix, sep="/")
  verbose = "1" if verbose else "0"
  libfn = get_lib_fn("antsRegistration")
  # Build arguments, see https://github.com/ANTsX/ANTsPy/blob/78dad33ca4e12ae605e80d8d990ec73a52abe273/ants/registration/interface.py
  args = [
    "--dimensionality", "3",
    "--float", "1",
    "-o", "[%s,%s,%s]" % (outprefix, outprefix + "Warped.nii.gz", outprefix + "InverseWarped.nii.gz"),
    "--interpolation", "Linear",
    "--use-histogram-matching", "1",
    
    # Rigid Transformation
    #     --transform Rigid[0.1] \
    #     --metric MI[fixedImage, movingImage, 1, 32, Random, 0.25] \
    #     --convergence [1000x500x250x0,1e-7,10] \
    #     --shrink-factors 12x8x4x1 \
    #     --smoothing-sigmas 5x4x3x1vox \
    "--transform", "Rigid[2.0]",
    "--metric", "MI[%s,%s,1,32,Random,0.25]" % (fixed_path, moving_path),
    "--convergence", "[1000x500x250x0,1e-7,10]",
    "--shrink-factors", "12x8x4x1",
    "--smoothing-sigmas", "5x4x3x1vox",

    # Affine Transformation
    #   --transform Affine[0.1] \
    #   --metric MI[/path/to/MNI152b.nii.gz,/path/to/subject_T1w.nii.gz,1,32,Random,0.25] \
    #   --convergence [1000x500x250x0,1e-7,10] \
    #   --shrink-factors 12x8x4x1 \
    #   --smoothing-sigmas 5x4x3x1vox \
    "--transform", "Affine[0.1]",
    "--metric", "MI[%s,%s,1,32,Random,0.25]" % (fixed_path, moving_path),
    "--convergence", "[1000x500x250x0,1e-7,10]",
    "--shrink-factors", "12x8x4x1",
    "--smoothing-sigmas", "5x4x3x1vox",

    # SyN Transformation
    #   --transform SyN[0.3,4,3] \
    #   --metric MI[/path/to/MNI152b.nii.gz,/path/to/subject_T1w.nii.gz,1,32,Random,0.25] \
    #   --convergence [1000x500x500x0,1e-7,7] \
    #   --shrink-factors 8x4x4x1 \
    #   --smoothing-sigmas 4x3x1x1vox
    "--transform", "SyN[0.3,4.0,3.0]",
    "--metric", "MI[%s,%s,1,32,Random,0.25]" % (fixed_path, moving_path),
    "--convergence", "[1000x500x500x0,1e-7,7]",
    "--shrink-factors", "8x4x4x1",
    "--smoothing-sigmas", "4x3x1x1vox",
    
    # subcortical refine stage
    "--transform", "SyN[0.3,4,3]",
    "--metric", "MI[%s,%s,1,32,Random,0.25]" % (fixed_path, moving_path),
    "--convergence", "[200x50x10x0,1e-6,7]",
    "--shrink-factors", "4x4x2x1",
    "--smoothing-sigmas", "2x2x1x1vox",

    "-v", verbose
  ]
  # Clean args
  processed_args = ants_process_arguments(args)
  reg_exit = libfn(processed_args)
  if verbose == "1":
    print("---- Call arguments --------------------------------")
    print("antsRegistration %s" % " ".join(processed_args))
  if not reg_exit == 0:
    raise RuntimeError(f"Registration failed with error code {reg_exit}")
  return {
    "outprefix"    : outprefix,
    "warpedmovout" : outprefix + "Warped.nii.gz",
    "warpedfixout" : outprefix + "InverseWarped.nii.gz",
    "transforms"   : [outprefix + "0GenericAffine.mat"]
  }

