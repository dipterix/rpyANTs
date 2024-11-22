
def get_lib_fn(fn):
  try:
    from ants import utils
    libfn = utils.get_lib_fn(fn)
  except:
    from ants.internal import get_lib_fn
    libfn = get_lib_fn(fn)
  return libfn


def ants_process_arguments(args):
  try:
    from ants import utils
    ret = utils._int_antsProcessArguments(args)
  except:
    from ants.internal import process_arguments
    ret = process_arguments(args)
  return ret