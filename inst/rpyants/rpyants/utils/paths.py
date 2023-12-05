import os
import re

def normalize_path(path, sep="/"):
  path = os.path.normpath(path)
  return sep.join(re.split(r"[\\/]+", path))
