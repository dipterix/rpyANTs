import ants
from typing import Union
from ..utils.paths import normalize_path, file_path
from ..utils.cache import _, StageContext, stage_image

def as_ANTsImage(x: Union[str, ants.core.ants_image.ANTsImage]):
    if isinstance(x, ants.core.ants_image.ANTsImage):
        return x
    if isinstance(x, str):
        img_path = normalize_path(x, sep="/")
        out = ants.image_read(img_path)
        return out
    raise TypeError("as_ANTsImage: `x` must be either a string or ANTsImage")

def brainmask(
        img_path: Union[str, ants.core.ants_image.ANTsImage],
        template_path: Union[str, ants.core.ants_image.ANTsImage],
        fixing_mask_path: Union[str, ants.core.ants_image.ANTsImage],
        working_path: str,
        verbose: bool = True):
    """
    Skullstrip an image using ANTsPyx

    Args:
        img (str): Path to the image to skullstrip

    Returns:
        str: Path to the skullstripped image
    """
    input_image=as_ANTsImage(img_path)
    fixing_img=as_ANTsImage(template_path)
    fixing_mask=as_ANTsImage(fixing_mask_path)
    
    # Need to get brain mask and find initial registration to template
    # Register brain using `SyNabp`: SyN optimized for abpBrainExtraction.
    with StageContext("SyNabp", "registration", working_path) as (ctx, res_abp):
        if res_abp is None:
            res_abp = ants.registration(
                fixed=fixing_img, moving=input_image,
                type_of_transform='SyNabp', verbose=verbose)
            ctx.result = res_abp
    # extract brain 
    with StageContext("brainmask", "image", working_path) as (ctx, brain_mask):
        if brain_mask is None:
            brain_mask = ants.apply_transforms(
                fixed=input_image, 
                moving=fixing_mask, 
                transformlist=res_abp["invtransforms"], 
                interpolator="nearestNeighbor")
            ctx.result = brain_mask
    return brain_mask
