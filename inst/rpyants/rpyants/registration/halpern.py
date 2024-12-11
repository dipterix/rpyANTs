# fixed_path = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/CT_RAW.nii.gz'
# moving_path = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/MRI_RAW.nii.gz'
# outprefix = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/test_'
# fixed = ants.image_read(fixed_path)
# moving = ants.image_read(moving_path)

import ants
from ..utils.paths import normalize_path
from ..utils.internals import get_lib_fn, ants_process_arguments

def halpern_coregister_ct_mri(fixed_path, moving_path, outprefix, verbose=True):
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
  # The difference is we fit coreg twice
  args = [
    "--dimensionality", "3",
    "--float", "0",
    "-o", "[%s,%s]" % (outprefix, outprefix + "Warped.nii.gz"),
    "--interpolation", "Linear",
    # "--winsorize-image-intensities", "[0.005, 0.995]",
    "--use-histogram-matching", "0",
    "--initial-moving-transform", "[%s,%s,1]" % (fixed_path, moving_path),
    "--transform", "Rigid[0.1]",
    "--shrink-factors", "8x4x2x1",
    "--metric", "MI[%s,%s,1,32,Regular,0.25]" % (fixed_path, moving_path),
    "--convergence", "[1000x500x250x100,1e-6,10]",
    "--shrink-factors", "8x4x2x1",
    "--smoothing-sigmas", "3x2x1x0vox",
    "--transform", "Rigid[0.1]",
    "--metric", "MI[%s,%s,1,32,Regular,0.25]" % (fixed_path, moving_path),
    "--convergence", "[1000x500x250x100,1e-6,10]",
    "--smoothing-sigmas", "3x2x1x0vox",
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
    "transforms"   : [outprefix + "0GenericAffine.mat"]
  }


# def halpern_coregister_template_mri(fixed, moving, outprefix, verbose=True):
#   '''
#   Register MRI to template for Casey Halpern's lab.
#   '''
#   # Make sure images are valid
#   if not isinstance(fixed, ants.core.ants_image.ANTsImage):
#     raise Exception("`fixed` is not an instance of `ANTsImage`")
#   if not isinstance(moving, ants.core.ants_image.ANTsImage):
#     raise Exception("`moving` is not an instance of `ANTsImage`")
#   # use forward slash to avoid escaping issues on windows
#   outprefix = normalize_path(outprefix, sep="/")
#   verbose = "1" if verbose else "0"
#   from ants import utils
#   fixed = fixed.clone("float")
#   moving = moving.clone("float")
#   # Do not use memory pointers
#   # # get memory pointers
#   # f = utils.get_pointer_string(fixed)
#   # m = utils.get_pointer_string(moving)
#   libfn = utils.get_lib_fn("antsRegistration")
#   # initial affine register
#   AffBasename = 't1_to_MNI_aff_'
#   ants.
# ANTs_command = AntsPth + str('antsRegistrationSyN.sh -d 3 -m %s -f %s -o %s -n %s -t a -r 2 -j 1') % (fnBrainNii, fnMNI, OutputAff, numcores)
# #print(ANTs_command)

# #os.system("sbatch -p normal -t 24:00:00 --mem=8192 --wrap=\"" + ANTs_command + "\"")

# SyNbasename = 't1_to_MNI_opt_';
# AffBrainNii = str(dpReg + str(str(AffBasename) + 'Warped.nii.gz'))
# OutNonl = str(dpReg + SyNbasename)
# cmd = AntsPth + str('antsRegistrationSyN.sh -d 3 -m %s -f %s -o %s -n %s -t b -r 2 -j 1 -x %s') % (AffBrainNii, fnMNI, OutNonl, numcores, fnVDCmask)

# #os.system("sbatch -p normal -t 24:00:00 --mem=8192 --wrap=\"" + cmd + "\"")
#   args = [
#     "--dimensionality", "3",
#     "--float", "0",
#     "-o", "[%s,%s]" % (outprefix, outprefix + "Warped.nii.gz"),
#     "--interpolation", "Linear",
#     # "--winsorize-image-intensities", "[0.005, 0.995]",
#     "--use-histogram-matching", "0",
#     "--initial-moving-transform", "[%s,%s,1]" % (fixed_path, moving_path),
#     "--transform", "Rigid[0.1]",
#     "--shrink-factors", "8x4x2x1",
#     "--metric", "MI[%s,%s,1,32,Regular,0.25]" % (fixed_path, moving_path),
#     "--convergence", "[1000x500x250x100,1e-6,10]",
#     "--shrink-factors", "8x4x2x1",
#     "--smoothing-sigmas", "3x2x1x0vox",
#     "--transform", "Rigid[0.1]",
#     "--metric", "MI[%s,%s,1,32,Regular,0.25]" % (fixed_path, moving_path),
#     "--convergence", "[1000x500x250x100,1e-6,10]",
#     "--smoothing-sigmas", "3x2x1x0vox",
#     "-v", verbose
#   ]
#   # Clean args
#   processed_args = utils._int_antsProcessArguments(args)
#   reg_exit = libfn(processed_args)
#   if verbose == "1":
#     print("---- Call arguments --------------------------------")
#     print("antsRegistration %s" % " ".join(processed_args))
#   if not reg_exit == 0:
#     raise RuntimeError(f"Registration failed with error code {reg_exit}")
#   return {
#     "outprefix"    : outprefix,
#     "warpedmovout" : outprefix + "Warped.nii.gz",
#     "transforms"   : [outprefix + "0GenericAffine.mat"]
#   }
