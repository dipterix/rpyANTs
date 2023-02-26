__all__ = [
    "get_pointer_string",
    "_int_antsProcessArguments",
]

import ants.utils as utils
from ants import lib
from ants.core import ants_image as iio
import os

def get_tmp_nii_file():
  try:
    # Use rpyANTs injected method
    tmpfile = utils.rpyANTsInjection.temp_nii_file()
  except:
    # standalone
    import tempfile
    tmpfile = tempfile.mktemp(suffix=".nii.gz")
  return os.path.normpath(tmpfile).replace(os.sep, '/')

def get_pointer_string(image):
  tmpfile = get_tmp_nii_file()
  image.to_file(tmpfile)
  return tmpfile

def _int_antsProcessArguments(args):
  """
  Needs to be better validated.
  """
  p_args = []
  if isinstance(args, dict):
    for argname, argval in args.items():
      if "-MULTINAME-" in argname:
        # have this little hack because python doesnt support
        # multiple dict entries w/ the same key like R lists
        argname = argname[: argname.find("-MULTINAME-")]
      if argval is not None:
        if len(argname) > 1:
          p_args.append("--%s" % argname)
        else:
          p_args.append("-%s" % argname)
        if isinstance(argval, iio.ANTsImage):
          p_args.append(get_pointer_string(argval))
        elif isinstance(argval, list):
          for av in argval:
            if isinstance(av, iio.ANTsImage):
              av = get_pointer_string(av)
            elif str(arg) == "True":
              av = str(1)
            elif str(arg) == "False":
              av = str(0)
            p_args.append(av)
        else:
          p_args.append(str(argval))
  elif isinstance(args, list):
    for arg in args:
      if isinstance(arg, iio.ANTsImage):
        pointer_string = get_pointer_string(arg)
        p_arg = pointer_string
      elif arg is None:
        pass
      elif str(arg) == "True":
        p_arg = str(1)
      elif str(arg) == "False":
        p_arg = str(0)
      else:
        p_arg = str(arg)
      p_args.append(p_arg)
  return p_args


setattr(utils, "get_pointer_string", get_pointer_string)
setattr(utils, "_int_antsProcessArguments", _int_antsProcessArguments)
setattr(utils, "rpyANTsInjected", True)
