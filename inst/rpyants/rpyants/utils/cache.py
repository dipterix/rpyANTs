import os
import ants
from .paths import file_path, ensure_basename, file_copy, file_move, normalize_path, unlink

class DummyVariable():
  def __init__(self):
    super().__setattr__("_", None)
  def __setattr__(self, name, value):
    pass

_ = DummyVariable()

def stage_image(img, name, root):
  path = file_path(root, name)
  ensure_basename(path)
  img.to_file(path)

def restore_image(name, root, strict = True):
  path = file_path(root, name)
  if not os.path.exists(path):
    if strict:
      raise Exception(f"Unable to restore image `{ name }`")
    else:
      return None
  return ants.image_read( path )

class StageContext():
  def __init__(self, prefix, result_type, root, verbose = True):
    self.prefix = prefix
    self.result_type = result_type
    self.root = root
    self.result = None
    self.needs_store = True
    self.verbose = verbose
    if result_type.startswith("image["):
      n_images = int(result_type[6:-1])
      self._n_images = n_images
  
  def _retrieve_result(self):
    root = self.root
    prefix = self.prefix
    result_type = self.result_type
    result = None
    if result_type == "registration":
      abs_prefix = file_path(root, prefix)
      result = {
        'warpedmovout' : restore_image(f"{ prefix }_warped_mov.nii.gz", root = root, strict = True),
        'fwdtransforms': [],
        'invtransforms': []
      }
      abs_prefix = file_path(root, prefix)
      # not gonna be 10...
      for ii in range(10):
        fwd_transform = f"{ abs_prefix }_fwdtransforms_ants{ ii }.nii.gz"
        inv_transform = f"{ abs_prefix }_invtransforms_ants{ ii }.nii.gz"
        if not os.path.exists(fwd_transform):
          fwd_transform = f"{ abs_prefix }_fwdtransforms_ants{ ii }.mat"
        if not os.path.exists(inv_transform):
          inv_transform = f"{ abs_prefix }_invtransforms_ants{ ii }.mat"
        if os.path.exists(fwd_transform) and os.path.exists(inv_transform):
          result['fwdtransforms'].append(fwd_transform)
          result['invtransforms'].append(inv_transform)
        else:
          break
      if len(result['fwdtransforms']) == 0:
        raise Exception("No transform data found")
    elif result_type.startswith("image["):
      # suffix must be from 0 to n_images-1, starting from 0
      result = [restore_image(name=f"{ prefix }_{ ii }.nii.gz", root = root, strict = True) for ii in range(self._n_images)]
    else:
      result = restore_image(name=f"{ prefix }.nii.gz", root = root, strict = True)
    if result is not None:
      self.result = result
    return result
  
  def _store_result(self, result):
    if result is None:
      return 
    root = self.root
    prefix = self.prefix
    result_type = self.result_type
    unlink_list = []
    if result_type.startswith("image["):
      if self._n_images <= 0:
        return
      if not isinstance(result, (tuple, list)):
        return
    if result_type == "registration":
      abs_prefix = file_path(root, prefix)
      warpedmovout = result['warpedmovout']
      if isinstance(warpedmovout, ants.core.ants_image.ANTsImage):
        stage_image(img=result['warpedmovout'], name=f"{ prefix }_warped_mov.nii.gz", root=root)
      else:
        warpedmovout = normalize_path(warpedmovout)
        unlink_list.append(warpedmovout)
        to_path = f"{ abs_prefix }_warped_mov.nii.gz"
        file_copy(warpedmovout, to_path)
        result['warpedmovout'] = to_path
      # stage_image(result['warpedmovout'], f"{ prefix }_warped_mov.nii.gz")
      for name in ["fwdtransforms", "invtransforms"]:
        for ii in range(len(result[name])):
          from_path = result[name][ii]
          unlink_list.append(from_path)
          if from_path.lower().endswith("nii.gz"):
            to_path = f"{ abs_prefix }_{ name }_ants{ ii }.nii.gz"
          else:
            to_path = f"{ abs_prefix }_{ name }_ants{ ii }.mat"
          file_copy(from_path, to_path)
          result[name][ii] = to_path
    elif result_type.startswith("image["):
      # suffix must be from 1 to n_images, not starting from 0
      for ii in range( self._n_images ):
        stage_image(img=result[ ii ], name=f"{ prefix }_{ ii }.nii.gz", root = root)
    else:
      stage_image(img=result, name=f"{ prefix }.nii.gz", root = root)
    for path in unlink_list:
      if os.path.exists(path):
        unlink(path)
  
  def __enter__(self):
    self.result = None
    try:
      self.result = self._retrieve_result()
      if self.result is not None:
        self.needs_store = False
        if self.verbose:
          print(f"Cache found for data `{ self.prefix }`")
    except Exception as e:
      if self.verbose:
        print(f"Cache not found for data `{ self.prefix }` (reason: { e })")
    return (self, self.result)
  
  def __exit__(self, *args, **kwargs):
    if self.needs_store and self.result is not None:
      if self.verbose:
        print(f"Updating data `{ self.prefix }`")
      self._store_result(self.result)
    self.result = None

