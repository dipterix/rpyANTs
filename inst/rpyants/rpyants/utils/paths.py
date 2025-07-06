import os
import re
import shutil
from typing import Union

def r_user_data_dir() -> str:
  package = "rpyANTs"
  which = "data"
  re = os.environ.get("R_USER_DATA_DIR", None)
  if re is None:
    re = os.environ.get("XDG_DATA_HOME", None)
  if re is None:
    os_type = "unknown"
    try:
      import platform
      os_type = platform.system().lower()
    except ImportError:
      if os.name == "nt":
        os_type = "windows"
      elif os.name == "posix":
        # check if file_path(home, "Library", "Application Support", "org.R-project.R") exists
        mac_data_path = os.path.join(os.path.expanduser("~"), "Library", "Application Support", "org.R-project.R")
        if os.path.exists(mac_data_path):
          os_type = "darwin"
        else:
          os_type = "linux"
    
    if os_type == "windows":
      appdata = os.environ.get("APPDATA", None)
      if appdata is not None:
        re = os.path.join(appdata, "R", "data")
    elif os_type == "darwin":
      re = os.path.join(os.path.expanduser("~"), "Library", "Application Support", "org.R-project.R")
    if re is None:
      re = os.path.join(os.path.expanduser("~"), ".local", "share")
  re = os.path.join(re, "R", package)
  return re

def antspynet_cache_dir() -> str:
  return os.path.join(r_user_data_dir(), "keras", "ANTsXNet")

def try_import_antspynet():
  try:
    import antspynet
    cache_dir = ensure_dir(antspynet_cache_dir())
    if os.path.exists(cache_dir):
      antspynet.set_antsxnet_cache_directory(cache_dir)
    return antspynet
  except ImportError:
    return None

def normalize_path(path, sep="/"):
  path = os.path.normpath(path)
  return sep.join(re.split(r"[\\/]+", path))

def unlink(path: str) -> bool:
  '''
  Remove a file or directory.

  @param path: The path to the file or directory.
  @type path: str

  @return: True if the file or directory was removed, False otherwise.
  @rtype: bool
  '''
  if os.path.exists(path):
    if os.path.isdir(path):
      shutil.rmtree(path, ignore_errors=True)
    else:
      os.unlink(path)
    return True
  return False

def ensure_dir(path: str) -> str:
  '''
  Ensure that the directory exists. If not, create it.

  @param path: The path to the directory.
  @type path: str

  @return: The path to the directory.
  @rtype: str
  '''
  if not os.path.exists(path):
    os.makedirs(path)
  return normalize_path(path)

def ensure_basename(filepath: Union[str, None]) -> str:
  '''
  Ensure that the directory of the file exists. If not, create it.

  @param filepath: The path to the file.
  @type filepath: str

  @return: The path to the file.
  @rtype: str
  '''
  if filepath is None:
    raise TypeError("Unable to find basename of path `None`")
  ensure_dir(os.path.dirname(filepath))
  return normalize_path(filepath)

def file_path(*args) -> str:
  '''
  Return the normalized path to a file.

  @param *args: The path components.
  @type *args: str

  @return: The normalized path to the file.
  @rtype: str
  '''
  return normalize_path(os.path.join(*args))

def parse_bids_filename(path):
  '''
  Parse a BIDS filename.

  @param path: The path to the BIDS filename.
  @type path: str

  @return: A dictionary containing the parsed BIDS filename.
  @rtype: dict
  '''
  # Remove the path and extension
  filename = os.path.basename(path)
  # Split the filename into components
  components = filename.split("_")
  fname = components.pop(-1)
  file_type_comp = fname.split(".")
  file_type = file_type_comp.pop(0)
  file_ext = ".".join(file_type_comp)
  parsed = {
    'fullname': filename,
    'name': fname,
    'ext': file_ext.lower(),
    'type': file_type,
    'components': {},
  }
  last_entity = None
  for i, c in enumerate(components):
    if c.find("-") == -1:
      k, v = last_entity, c
    else:
      k, v = c.split("-")
      last_entity = k
    if k is not None:
      if parsed['components'].get(k, None):
        old_v = parsed['components'][k]
        v = f"{ old_v }_{ v }"
      parsed['components'][k] = v
  return parsed


def file_copy(src, dst, auto_mkdirs : bool = True):
  '''
  Copy a file from `src` to `dst`.

  @param src: The source file.
  @type src: str

  @param dst: The destination file.
  @type dst: str

  @param auto_mkdirs: Automatically create the parent directories if they do not exist.
  @type auto_mkdirs: bool

  @return: The destination file.
  @rtype: str
  '''
  if auto_mkdirs:
    ensure_dir(os.path.dirname(dst))
  shutil.copyfile(src, dst)
  return normalize_path(dst)

def file_move(src, dst, auto_mkdirs : bool = True):
  '''
  Move a file from `src` to `dst`.

  @param src: The source file.
  @type src: str

  @param dst: The destination file.
  @type dst: str

  @param auto_mkdirs: Automatically create the parent directories if they do not exist.
  @type auto_mkdirs: bool

  @return: The destination file.
  @rtype: str
  '''
  if auto_mkdirs:
    ensure_dir(os.path.dirname(dst))
  shutil.move(src, dst)
  return normalize_path(dst)

def to_bids_prefix(components : dict, **kwargs) -> str:
  '''
  Convert a parsed BIDS filename to a BIDS file name prefix.
  For example, { 'sub': 'AnonSEEG', 'ses': 'postop', 'space': 'MNI152NLin2009bAsym', 'rec': 'tonemapped', 'desc': 'preproc' }
  will return `sub-AnonSEEG_ses-postop_space-MNI152NLin2009bAsym_rec-tonemapped_desc-preproc_`

  @param components: The parsed BIDS components (see `components` in `parse_bids_filename`).
  @type parsed: dict

  @param kwargs: Additional components to add to the BIDS filename, will overwrite existing components.
  @type kwargs: dict

  @return: The BIDS filename.
  @rtype: str
  '''
  parts = []
  reserved_keys = ['sub', 'ses', 'acq', 'rec', 'run', 'task', 'dir', 'space', 'from', 'to', 'desc']
  for k in reserved_keys:
    if k in kwargs:
      parts.append(f"{k}-{kwargs[k]}")
    elif k in components:
      parts.append(f"{k}-{components[k]}")
  for k, v in kwargs.items():
    if k not in reserved_keys:
      parts.append(f"{k}-{v}")
      reserved_keys.append(k)
  for k, v in components.items():
    if k not in reserved_keys:
      parts.append(f"{k}-{v}")
  return "_".join(parts) + "_"
