# ct_path = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/CT_RAW.nii.gz'
# mr_path = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/MRI_RAW.nii.gz'
# template_path = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/CT_Template.nii.gz'
# outprefix = '/Users/dipterix/rave_data/raw_dir/testtest/rave-imaging/coregistration/test_'
# fixed = ants.image_read(fixed_path)
# moving = ants.image_read(moving_path)

# import rpyants
# yael = rpyants.registration.YAELPreprocess("PAV038", "/Users/dipterix/rave_data/raw_dir/PAV038/rave-imaging")
# yael.get_native_mapping("CT")
import os
import json
from typing import Union
import numpy as np
import ants

from ..utils.paths import normalize_path, ensure_dir, file_path, parse_bids_filename, to_bids_prefix, file_copy, ensure_basename, unlink
from ..utils.transforms import ants_AffineTransform_to_m44, M44_LPS_to_RAS, get_xform, set_xform
from .halpern import halpern_coregister_ct_mri
from .normalization import normalize_to_template_syn, normalization_with_atropos



# There are N types of registrations:
# 1. CT, T2, ... to T1 MRI (rigid-body)
# 2. T1 MRI to Template (non-linear)
class YAELPreprocess():
    @property
    def work_path(self):
        return self._work_path
    
    def __init__(self, subject_code : str, work_path : str, image_types : Union[list, tuple] = ["CT", "T1w", "T2w", "FLAIR", "preopCT", "T1wContrast", "fGATIR", "DWI"] ):
        '''
        Initialize a YAELPreprocess object.
        
        @param subject_code: The subject code (with or without `sub-` prefix), for example, `sub-DemoSubject` or `DemoSubject`.
        @type subject_code: str

        @param work_path: The path to the working directory, for example, `rave-imaging/`.
        @type work_path: str

        @param image_types: Image types allowed to query
        @type image_types: list or tuple

        '''
        
        # List of allowed image types
        # CT: postop CT
        # T1w: preop T1w MRI without contrast
        # T2w: preop T2w MRI 
        # FLAIR: preop FLAIR MRI
        # preopCT: preop CT (for showing blood vessels)
        # T1wContrast: preop T1w MRI with contrast (for showing blood vessels)
        # fGATIR: preop fast Gray Matter Acquisition T1 Inversion Recovery
        self.allowed_image_types = [ *image_types ]
        # Normalization weights relative to T1w (weight=1)
        self.normalization_weights = {
            "CT"      : 0.60,
            "preopCT" : 0.80,
            "postopCT": 0.60,
            "FLAIR"   : 0.25,
            "postopFLAIR" : 0.2,
            "preopFLAIR"  : 0.3,
        }
        # Each item is a dictionary returned by `parse_bids_filename`
        self._images = {}

        self._subject_code = subject_code
        if subject_code.startswith("sub-"):
            self._subject_code = subject_code[4:]
        self._work_path = normalize_path( work_path )
        self.reload_images()
    
    def _fix_image_path(self, parsed : dict, name = None) -> Union[dict, None]:
        if parsed is None:
            return None
        if parsed.get('components', None) is None:
            return None
        parsed['components']['sub'] = self._subject_code
        parsed['components']['desc'] = 'preproc'
        if parsed['ext'] == "nii":
            parsed['ext'] = "nii.gz"
        if parsed['components'].get('ses', None) is None:
            image_type = parsed['type']
            if image_type.startswith("preop"):
                parsed['components']['ses'] = 'preop'
            elif image_type.startswith("preop"):
                parsed['components']['ses'] = 'postop'
            elif image_type in self.allowed_image_types:
                if parsed['type'] == "CT":
                    parsed['components']['ses'] = 'postop'
                else:
                    parsed['components']['ses'] = 'preop'
        if name is not None:
            parsed['name'] = name + "." + parsed['ext']
        else:
            parsed['name'] = parsed['type'] + "." + parsed['ext']
        parsed['fullname'] = to_bids_prefix(parsed['components']) + parsed['name']
        return parsed
    
    def reload_images(self): 
        paths = {}
        input_path = file_path(self._work_path, "inputs", "anat")
        if os.path.exists(input_path) and os.path.isdir(input_path):
            for filename in os.listdir(input_path):
                if filename.startswith("."):
                    continue
                if not os.path.isdir(file_path(input_path, filename)):
                    try:
                        parsed = parse_bids_filename(filename)
                        if parsed is not None and parsed['type'] in self.allowed_image_types:
                            parsed['folder'] = "inputs/anat"
                            paths[parsed['type']] = parsed
                    except:
                        pass
        self._images = paths
    
    def input_image_path(self, type : str, relative : bool = False) -> Union[str, None]:
        '''
        Get the path to the image.

        @param type: The type of the image (e.g., `CT`, `T1w`, "T2w", "FLAIR", "preopCT", "T1wContrast", "fGATIR")
        @type type: str

        @param relative: If True, return the relative path (to the working directory). 
                Otherwise return the absolute path.
        @type relative: bool

        @return: The path to the image.
        @rtype: str if found or None is missing
        '''
        item = self._images.get(type, None)
        if item is None:
            return None
        if relative:
            return file_path(item['folder'], item['fullname'])
        return file_path(self._work_path, item['folder'], item['fullname'])
    
    def format_path(self, folder : str, name : str, ext : str, relative : bool = False, **kwargs):
        '''
        Format the path to the file.

        @param folder: The folder where the file is expected to be found.
        @type folder: str

        @param relative: If True, return the relative path (to the working directory). 
                Otherwise return the absolute path.
        @type relative: bool

        @param kwargs: Additional arguments
        @type kwargs: dict
            * Other argument keys are BIDS components (e.g., `sub`, `ses`, `space`, `rec`, `desc`)

        @return: The path to the file.
        @rtype: str if found or None is missing
        '''
        parsed = {
            'folder': folder,
            'type' : name,
            'ext' : ext,
            'components': {},
        }
        parsed['components'].update(kwargs)
        parsed = self._fix_image_path(parsed)
        if parsed is None:
            return None
        if relative:
            return file_path(parsed['folder'], parsed['fullname'])
        return file_path(self._work_path, parsed['folder'], parsed['fullname'])
    def expected_image_path(self, type : str, folder : Union[str, None] = None, name = None, **kwargs):
        '''
        Get the expected path to the image.

        @param type: The type of the image (e.g., `CT`, `T1w`, "T2w", "FLAIR", "preopCT", "T1wContrast", "fGATIR")
        @type type: str

        @param folder: The folder where the image is expected to be found.
        @type folder: str

        @param kwargs: Additional arguments
        @type kwargs: dict
            * Other argument keys are BIDS components (e.g., `sub`, `ses`, `space`, `rec`, `desc`)

        @return: The path to the image.
        @rtype: str
        '''
        item = self._images.get(type, None)
        if item is None:
            return None
        parsed = {
            'folder': folder,
            'ext': item['ext'],
            'type': type,
            'components': {},
        }
        parsed['components'].update(kwargs)
        parsed = self._fix_image_path(parsed, name = name)
        if parsed is None:
            raise ValueError(f"Unable to get image path for `{name}`")
        return file_path(self._work_path, parsed['folder'], parsed['fullname'])
    
    def set_image(self, type : str, path : str, **kwargs):
        '''
        Set the path to the image.

        @param type: The type of the image (e.g., `CT`, `T1w`, "T2w", "FLAIR", "preopCT", "T1wContrast", "fGATIR")
        @type type: str

        @param path: The path to the image.
        @type path: str

        @param kwargs: Additional arguments
        @type kwargs: dict
            * overwrite : If True, overwrite the existing image. Otherwise throw an error if image exists.
            * Other argument keys are BIDS components (e.g., `sub`, `ses`, `space`, `rec`, `desc`)
        '''
        if not os.path.exists(path):
            raise FileNotFoundError(f"Invalid image path: {path}")
        if type not in self.allowed_image_types:
            # raise ValueError(f"Invalid image type: {type}. Must be one of {self.allowed_image_types}")
            self.allowed_image_types.append(type)
        overwrite = kwargs.get('overwrite', False)
        existing_path = self.input_image_path(type)
        if not overwrite and existing_path is not None:
            raise FileExistsError(f"Image type {type} already exists. If you want to overwrite it, set `overwrite=True`.")
        try:
            parsed = parse_bids_filename(os.path.basename(path))
            if parsed is None:
                raise ValueError("Invalid BIDS filename")
            parsed['folder'] = "inputs/anat"
            parsed['type'] = type
            parsed['ext'] = "nii.gz"
            parsed = self._fix_image_path(parsed)
        except:
            parsed = self._fix_image_path({
                'folder': "inputs/anat",
                'ext': "nii.gz",
                'type': type,
                'components': {
                    'sub': self._subject_code,
                },
            })
        self._images[type] = parsed
        if existing_path is not None:
            unlink(existing_path)
        # file_copy(path, self.input_image_path(type))
        target_path = self.input_image_path(type)
        if target_path is None:
            raise ValueError(f"Unable to determine path for image type `{ type }`")
        ants.image_read(path).to_file(ensure_basename(target_path))
    
    @property
    def input_ct_path(self):
        return self.input_image_path("CT")
    
    @property
    def input_t1w_path(self):
        return self.input_image_path("T1w")
    
    @property
    def input_t2w_path(self):
        return self.input_image_path("T2w")
    
    @property
    def input_flair_path(self):
        return self.input_image_path("FLAIR")
    
    @property
    def input_preopCT_path(self):
        return self.input_image_path("preopCT")
    
    @property
    def input_fGATIR_path(self):
        return self.input_image_path("fGATIR")
    
    @property
    def input_t1wcontrast_path(self):
        return self.input_image_path("T1wContrast")
    
    def register_to_T1w(self, type : str, reverse : bool= False, verbose : bool=True):
        '''
        Register an image to the T1w image.

        @param type: The type of the image (e.g., `CT`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR")
        @type type: str

        @param reverse: If True, register the T1w image to the image (switch the fixed and moving images).
        @type reverse: bool

        @param verbose: If True, print verbose output.
        @type verbose: bool
        '''
        if type == "T1w":
            raise ValueError("Cannot register T1w to itself.")
        fix_path = self.input_t1w_path
        if fix_path is None:
            raise FileNotFoundError(f"Missing image for type `T1w`. Please set the image first via `preprocess.set_image`.")
        mov_path = self.input_image_path(type)
        if mov_path is None:
            raise FileNotFoundError(f"Missing image for type `{type}`. Please set the image first via `preprocess.set_image`.")
        outprefix = file_path(self._work_path, "tmp", f"coregister_{type}_with_T1w", "out_")
        outputdir = os.path.dirname(outprefix)
        if os.path.exists(outputdir):
            unlink(outputdir)
        ensure_dir(outputdir)
        if reverse:
            fix_path, mov_path = mov_path, fix_path
        rigid_reg = halpern_coregister_ct_mri(fix_path, mov_path, outprefix, verbose=verbose)
        # copy the transform to the expected location
        transform_srcpath = rigid_reg['transforms'][0]
        transform_from_t1_prefix = to_bids_prefix({ 'sub': self._subject_code, 'from': 'T1w', 'to': type, 'desc': 'rigid', })
        transform_to_t1_prefix = to_bids_prefix({ 'sub': self._subject_code, 'from': type, 'to': 'T1w', 'desc': 'rigid', })
        affine = ants.read_transform(transform_srcpath)
        if reverse:
            file_copy(transform_srcpath, ensure_basename( file_path( self._work_path, "coregistration/transformations", transform_from_t1_prefix + "ants.mat" ) ))
            ants.write_transform(affine.invert(), ensure_basename( file_path( self._work_path, "coregistration/transformations", transform_to_t1_prefix + "ants.mat" ) ))
        else:
            file_copy(transform_srcpath, ensure_basename( file_path( self._work_path, "coregistration/transformations", transform_to_t1_prefix + "ants.mat" ) ))
            ants.write_transform(affine.invert(), ensure_basename( file_path( self._work_path, "coregistration/transformations", transform_from_t1_prefix + "ants.mat" ) ))
        # calculate the CT RAS to T1 RAS transform
        # ANTs affine is from CT LPS to MR LPS
        affine = ants.read_transform(transform_srcpath)
        fix_lps_to_mov_lps = ants_AffineTransform_to_m44(affine)
        fix_ras_to_mov_ras = np.matmul(M44_LPS_to_RAS, np.matmul(fix_lps_to_mov_lps, M44_LPS_to_RAS))
        if reverse:
            ct_ras_to_mr_ras = fix_ras_to_mov_ras
            mr_ras_to_ct_ras = np.linalg.inv(fix_ras_to_mov_ras)
        else:
            ct_ras_to_mr_ras = np.linalg.inv(fix_ras_to_mov_ras)
            mr_ras_to_ct_ras = fix_ras_to_mov_ras
        np.savetxt( 
            file_path( self._work_path, "coregistration/transformations", transform_to_t1_prefix + "ras2ras.tsv" ),
            ct_ras_to_mr_ras, fmt = "%f", delimiter = "\t", newline="\n"
        )
        np.savetxt(
            file_path( self._work_path, "coregistration/transformations", transform_from_t1_prefix + "ras2ras.tsv" ),
            mr_ras_to_ct_ras, fmt = "%f", delimiter = "\t", newline="\n"
        )
        # copy the output image to the expected location
        ct_path = self.input_image_path(type)
        # ct_img = nib.load(ct_path)
        # current_sform = ct_img.get_sform()
        ct_img = ants.image_read(ct_path)
        current_sform = get_xform(ct_img)
        new_sform = np.matmul(ct_ras_to_mr_ras, current_sform)
        np.savetxt(
            file_path( self._work_path, "coregistration/transformations", transform_to_t1_prefix + "vox2ras.tsv" ),
            new_sform, fmt = "%f", delimiter = "\t", newline="\n"
        )
        # ct_img.set_sform(new_sform)
        # ct_img.set_qform(new_sform)
        ct_img = set_xform(ct_img, new_sform)
        ct_img.to_filename( ensure_basename( self.expected_image_path( type, "coregistration/anat", space = "scanner" ) ) )
        # save the mapping configurations
        mappings = self.get_native_mapping(type, relative = True)
        log_file = ensure_basename( self.format_path(folder="coregistration/log", name = "mapping", ext="json", space = "scanner", orig = type) )
        with open(log_file, "w") as f:
            json.dump(mappings, f)
        return mappings

    def get_native_mapping(self, type : str = "CT", relative : bool = False):
        '''
        Get the mappings for the native image.

        @param type: The type of the image modality (e.g., "CT", "T2w", "FLAIR", "preopCT", "T1wContrast", "fGATIR")
        @type type: str

        @param relative: If True, return the relative path (to the working directory).
        @type relative: bool

        @return: Dictionary containing the mappings from this modality to `T1w` image.
        @rtype: dict
        '''
        if type not in self.allowed_image_types:
            raise ValueError(f"Invalid image type: {type}. Must be one of {self.allowed_image_types}")
        tranform_rootdir = file_path(self._work_path, "coregistration/transformations")
        aligned_rootdir = file_path(self._work_path, "coregistration/anat")
        if relative:
            transform_prefix = "coregistration/transformations"
            aligned_prefix = "coregistration/anat"
        else:
            transform_prefix = tranform_rootdir
            aligned_prefix = aligned_rootdir
        result = {
            'subject': self._subject_code,
            'work_path': self._work_path,
            f'{type}_path': self.input_image_path(type, relative = relative),
            'T1w_path': self.input_image_path("T1w", relative = relative),
            'mappings': {},
        }
        mappings = result['mappings']
        if not os.path.exists(tranform_rootdir):
            ensure_dir(tranform_rootdir)
        for filename in os.listdir( tranform_rootdir ):
            if os.path.isfile( file_path(tranform_rootdir, filename) ):
                try:
                    parsed = parse_bids_filename(filename)
                    parsed_components = parsed['components']
                    if parsed_components.get('from', None) == type and parsed_components.get('to', None) == "T1w":
                        mappings[parsed['type']] = file_path(transform_prefix, filename)
                except:
                    continue
        if not os.path.exists(aligned_rootdir):
            ensure_dir(aligned_rootdir)
        for filename in os.listdir( aligned_rootdir ):
            if os.path.isfile( file_path(aligned_rootdir, filename) ):
                try:
                    parsed = parse_bids_filename(filename)
                    parsed_components = parsed['components']
                    mtype = parsed.get('type', None)
                    if mtype == type:
                        mappings[f"{type}_in_T1w"] = file_path(aligned_prefix, filename)
                except:
                    continue
        return result

    def map_to_template(self, 
                        template_name : str, template_path : str, template_mask_path = None,
                        native_type : str = "T1w", native_mask_path = None,
                        use_images : Union[str, list] = "all",
                        verbose : bool=True, **kwargs):
        '''
        Register an image to the template.

        @param template_name: The name of the template image (e.g., `MNI152NLin2009bAsym`)
        @type template_name: str

        @param template_path: The path to the template image.
        @type template: str

        @param template_mask_path: The path to the template brain mask
        @type template: str

        @param native_type: The type of the image (e.g., `CT`, `T1w`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR")
        @type native_type: str

        @param native_mask_path: The path to the native underlay brain mask
        @type template: str

        @param use_images: Images to use for normalization besides underlaying images. This can be "all" for all the image types, "none" for none, or a list of image types
        @type use_images: str or list of image types

        @param verbose: If True, print verbose output.
        @type verbose: bool
        
        @param **kwargs: additional parameters passed to `normalization_with_atropos`
        '''
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Invalid template path: {template_path}")
        if native_type not in self.allowed_image_types:
            raise ValueError(f"Invalid image type: {native_type}. Must be one of {self.allowed_image_types}")
        native_path = self.input_image_path(native_type)
        if native_path is None:
            raise FileNotFoundError(f"Missing image for type `{native_type}`. Please set the image first via `preprocess.set_image`.")
        # Check the image registrations for other modalities
        other_modalities = []
        other_images = {}
        if use_images is not None:
            if isinstance(use_images, str):
                if use_images == "all":
                    other_modalities = [*self.allowed_image_types]
                else:
                    other_modalities = [use_images]
            else:
                other_modalities = [*use_images]
        other_modalities = [x for x in other_modalities if x in self.allowed_image_types]
        for modality in other_modalities:
            if modality == native_type:
                continue
            mapping=self.get_native_mapping(type=modality, relative=False)
            mapped_path=mapping['mappings'].get(f"{modality}_in_{native_type}", None)
            # register to T1 if the mapping is missing
            if native_type == "T1w" and mapped_path is None:
                orig_path=mapping.get(f"{modality}_path")
                if orig_path is not None:
                    self.register_to_T1w(type=modality, verbose=verbose)
                    mapping=self.get_native_mapping(type=modality, relative=False)
                    mapped_path=mapping['mappings'].get(f"{modality}_in_{native_type}", None)
            if mapped_path is not None:
                other_images[modality] = mapped_path
        fix_path, mov_path = template_path, native_path
        ants_outputdir = file_path(self._work_path, "tmp", f"coregister_{native_type}_with_{template_name}")
        resultdir = file_path(self._work_path, "normalization")
        if os.path.exists(ants_outputdir):
            unlink(ants_outputdir)
        ensure_dir(ants_outputdir)
        ensure_dir(resultdir)
        if template_mask_path is not None and os.path.exists(template_mask_path):
            template_mask = ants.image_read(template_mask_path)
        else:
            template_mask = None
        if native_mask_path is not None and os.path.exists(native_mask_path):
            native_mask = ants.image_read(native_mask_path)
        else:
            native_mask = None
        # for saving the registration result
        transform_from_template_prefix = to_bids_prefix({ 'sub': self._subject_code, 'from': template_name, 'to': native_type, 'desc': 'affine+SyN', })
        transform_to_template_prefix = to_bids_prefix({ 'sub': self._subject_code, 'from': native_type, 'to': template_name, 'desc': 'affine+SyN', })
        mov_paths = [mov_path]
        weights = [1.25]
        for modality, path in other_images.items():
            mov_paths.append( path )
            weights.append( self.normalization_weights.get(modality, 0.9) )
        # if True:
        #     from rpyants.registration.normalization import normalize_to_template_syn
        #     from rpyants.utils.paths import file_path
        #     fix_path='/Users/dipterix/rave_data/raw_dir/PAV038/rave-imaging/normalization/anat/MNI152NLin2009bAsym_orig.nii.gz'
        #     mov_paths=['/Users/dipterix/rave_data/raw_dir/PAV038/rave-imaging/inputs/anat/sub-PAV038_ses-preop_desc-preproc_T1w.nii.gz']
        #     weights = [1.25]
        #     ants_outputdir='/Users/dipterix/rave_data/raw_dir/PAV038/rave-imaging/tmp/coregister_T1w_with_MNI152NLin2009bAsym'
        #     verbose=True
        deformable_reg = normalization_with_atropos(
            fix_path=fix_path,
            mov_paths=mov_paths,
            working_path = ants_outputdir,
            weights = weights,
            verbose = verbose, 
            **kwargs
        )
        # save forward (mov -> fix) transforms
        # simple experiment to check the composition of the transforms: the following two are the same
        # np.matmul(ants_AffineTransform_to_m44(affine_1), ants_AffineTransform_to_m44(affine_0))
        # ants_AffineTransform_to_m44(ants.compose_ants_transforms([affine_0, affine_1]))
        file_copy( deformable_reg['fwdtransforms'][1], ensure_basename(
            file_path( self._work_path, "normalization/transformations", transform_to_template_prefix + "ants0.mat" )
        ))
        file_copy( deformable_reg['fwdtransforms'][0], ensure_basename(
            file_path( self._work_path, "normalization/transformations", transform_to_template_prefix + "ants1.nii.gz" )
        ))
        # calculate and save the inverse transforms (fix -> mov)
        file_copy( deformable_reg['invtransforms'][1], ensure_basename(
            file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants0.nii.gz" )
        ))
        ants.write_transform(
            ants.read_transform(deformable_reg['fwdtransforms'][1]).invert(),
            file_path( self._work_path, "normalization/transformations", transform_from_template_prefix + "ants1.mat" )
        )
        # help(ants.apply_transforms)
        # The commented code reproduces the same result as `deformable_reg['warpedmovout']`
        # warped_T1_img = ants.apply_transforms(
        #     fixed = fix_img, moving = mov_img,
        #     transformlist=[deformable_reg['fwdtransforms'][0], deformable_reg['fwdtransforms'][1], affine_reg['fwdtransforms'][0]],
        #     whichtoinvert = [False, False, False]
        # )
        # warped_T1_img.to_file(ensure_basename(
        #     self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "deformable")
        # ))
        ###
        # deformable_reg['warpedmovout'].to_file(ensure_basename(
        #     self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "deformable")
        # ))
        if isinstance(deformable_reg['warpedmovout'], str):
            file_copy(
                src = deformable_reg['warpedmovout'], 
                dst = self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "deformable")
            )
        else:
            normalization_morph_path = self.expected_image_path(native_type, "normalization/anat", name = native_type, space = template_name, transform = "deformable")
            deformable_reg['warpedmovout'].to_file(ensure_basename(normalization_morph_path))
        ###
        # fix_img = ants.image_read(fix_path)
        # mov_img = ants.image_read(mov_path)
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
        if isinstance(deformable_reg['warpedfixout'], str):
            file_copy(
                src = deformable_reg['warpedfixout'], 
                dst = self.expected_image_path(native_type, "normalization/anat", name = template_name, space = native_type, transform = "deformable")
            )
        else:
            normalization_morph_path = self.expected_image_path(native_type, "normalization/anat", name = template_name, space = native_type, transform = "deformable")
            deformable_reg['warpedfixout'].to_file(ensure_basename(normalization_morph_path))
        mappings = self.get_template_mapping(template_name = template_name, native_type = native_type, relative = True)
        if mappings is not None:
            log_file = ensure_basename( self.format_path(folder="normalization/log", name="mappings", ext="json", 
                                                        template = template_name, native = native_type) )
            with open(log_file, "w") as f:
                json.dump(mappings, f)
        return mappings

    def get_template_mapping(self, template_name : str, native_type : str = "T1w", relative : bool = False):
        '''
        Get the mapping from the native image to the template.

        @param template_name: The name of the template image (e.g., `MNI152NLin2009bAsym`)
        @type template_name: str

        @param native_type: The type of the image (e.g., `CT`, `T1w`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR"), default is `T1w`.
        @type native_type: str

        @param relative: If True, return the relative path (to the working directory).
        @type relative: bool

        @return: Dictionary containing the mapping from/to the native image to/from the template.
        @rtype: dict
        '''
        if native_type not in self.allowed_image_types:
            raise ValueError(f"Invalid image type: {native_type}. Must be one of {self.allowed_image_types}")
        tranform_rootdir = file_path(self._work_path, "normalization/transformations")
        if not os.path.exists(tranform_rootdir):
            return None
        forward_list = []
        inverse_list = []
        for filename in os.listdir( tranform_rootdir ):
            if os.path.isfile( file_path(tranform_rootdir, filename) ):
                try:
                    parsed = parse_bids_filename(filename)
                    parsed_components = parsed['components']
                    order = parsed['type'].replace("ants", "")
                    order = int(order)
                except:
                    continue
                if parsed_components.get('sub', None) == self._subject_code:
                    if order is not None:
                        if parsed_components.get('from', None).lower() == native_type.lower() and parsed_components.get('to', None).lower() == template_name.lower():
                            forward_list.append( (order, filename) )
                        elif parsed_components.get('from', None).lower() == template_name.lower() and parsed_components.get('to', None).lower() == native_type.lower():
                            inverse_list.append( (order, filename) )
        if len( forward_list ) == 0 or len( inverse_list ) == 0:
            return None
        forward_list.sort(key = lambda x: x[0], reverse = True)
        inverse_list.sort(key = lambda x: x[0], reverse = True)
        if relative:
            path_prefix = "normalization/transformations"
        else:
            path_prefix = tranform_rootdir
        return {
            'type': native_type,
            'template': template_name,
            'using': 'ants.apply_transforms',
            'native_to_template': {
                'transformlist' : [ file_path(path_prefix, x[1]) for x in forward_list ],
                'whichtoinvert' : [ False for _ in forward_list ]
            },
            'template_to_native': {
                'transformlist' : [ file_path(path_prefix, x[1]) for x in inverse_list ],
                'whichtoinvert' : [ False for _ in inverse_list ]
            }
        }
    
    def transform_image_from_template(self, path : Union[str, ants.core.ants_image.ANTsImage], template_name : str, native_type : str = "T1w", interpolator="nearestNeighbor", verbose : bool=True):
        '''
        Map the image from the template to the native image.

        @param path: The path to the image.
        @type path: str | ANTsImage

        @param template_name: The name of the template image (e.g., `MNI152NLin2009bAsym`)
        @type template_name: str

        @param native_type: The type of the image (e.g., `CT`, `T1w`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR"), default is `T1w`.
        @type native_type: str

        @return: The mapped image.
        @rtype: ANTsImage
        '''
        if isinstance(path, ants.core.ants_image.ANTsImage):
            moving_img = path
        else:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Invalid image path: {path}")
            moving_img = ants.image_read(path)
        map_info = self.get_template_mapping(template_name = template_name, native_type = native_type)
        if map_info is None:
            raise FileNotFoundError(f"Missing mapping from {native_type} to {template_name}. Please register the image to the template first.")
        map_args = map_info['template_to_native']
        mapped_img = ants.apply_transforms(
            fixed = ants.image_read(self.input_image_path(native_type)),
            moving = moving_img,
            interpolator = interpolator,
            verbose = verbose,
            **map_args
        )
        return mapped_img
    
    def transform_image_to_template(self, path : Union[str, ants.core.ants_image.ANTsImage], template_name : str, template_path : Union[str, ants.core.ants_image.ANTsImage], native_type : str = "T1w", interpolator="nearestNeighbor", verbose : bool=True):
        '''
        Map the image from the native image to the template.

        @param path: The path to the native image.
        @type path: str | ANTsImage

        @param template_name: The name of the template image (e.g., `MNI152NLin2009bAsym`)
        @type template_name: str

        @param template_path: The path to the template image
        @type template_path: str | ANTsImage

        @param native_type: The type of the image (e.g., `CT`, `T1w`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR"), default is `T1w`.
        @type native_type: str

        @return: The mapped image.
        @rtype: ANTsImage
        '''
        map_info = self.get_template_mapping(template_name = template_name, native_type = native_type)
        if map_info is None:
            raise FileNotFoundError(f"Missing mapping from {native_type} to {template_name}. Please register the image to the template first.")
        if isinstance(path, ants.core.ants_image.ANTsImage):
            moving_img = path
        else:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Invalid image path: {path}")
            moving_img = ants.image_read(path)
        if isinstance(template_path, ants.core.ants_image.ANTsImage):
            template_img = template_path
        else:
            if not os.path.exists(template_path):
                raise FileNotFoundError(f"Invalid template path: {template_path}")
            template_img = ants.image_read(template_path)
        map_args = map_info['native_to_template']
        mapped_img = ants.apply_transforms(
            fixed = template_img,
            moving = moving_img,
            interpolator = interpolator,
            verbose = verbose,
            **map_args
        )
        return mapped_img
    
    def transform_points_to_template(self, points : np.ndarray, template_name : str, native_type : str = "T1w", verbose : bool=True):
        '''
        Map the points from the template to the native image.

        @param points: The points (`nx3`) in the native space (specified by `native_type`), `RAS` (right-anterior-superior) coordinate system.
        @type points: np.ndarray

        @param template_name: The name of the template image (e.g., `MNI152NLin2009bAsym`)
        @type template_name: str

        @param native_type: The type of the image (e.g., `CT`, `T1w`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR"), default is `T1w`.
        @type native_type: str

        @return: The mapped points.
        @rtype: np.ndarray
        @raise FileNotFoundError: If the mapping is missing.

        @example:
        ```
        import numpy as np
        points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        # preprocess.set_image("T1w", "/path/to/T1w.nii.gz")
        # preprocess.map_to_template(template_path = "/path/to/MNI152NLin2009bAsym.nii.gz", template_name = "MNI152NLin2009bAsym")
        mapped_points = preprocess.map_points_to_template(points, template_name = "MNI152NLin2009bAsym")
        ```

        '''
        map_info = self.get_template_mapping(template_name = template_name, native_type = native_type)
        if map_info is None:
            raise FileNotFoundError(f"Missing mapping from {native_type} to {template_name}. Please register the image to the template first.")
        map_args = map_info['template_to_native']
        import pandas as pd
        df = pd.DataFrame({
            'x' : - points[:, 0],
            'y' : - points[:, 1],
            'z' : points[:, 2]
        })
        mapped_points = ants.apply_transforms_to_points(
            dim = 3, points = df, 
            verbose = verbose,
            **map_args
        )
        return np.array([
            - mapped_points['x'].values,
            - mapped_points['y'].values,
            mapped_points['z'].values
        ]).transpose()

    def transform_points_from_template(self, points : np.ndarray, template_name : str, native_type : str = "T1w", verbose : bool=True):
        '''
        Map the points from the template to the native image.

        @param points: The points (`nx3`) in the template space (specified by `template_name`), `RAS` (right-anterior-superior) coordinate system.
        @type points: np.ndarray

        @param template_name: The name of the template image (e.g., `MNI152NLin2009bAsym`)
        @type template_name: str

        @param native_type: The type of the image (e.g., `CT`, `T1w`, `T2w`, "FLAIR", "preopCT", "T1wContrast", "fGATIR"), default is `T1w`.
        @type native_type: str

        @return: The mapped points.
        @rtype: np.ndarray
        @raise FileNotFoundError: If the mapping is missing.

        @example:
        ```
        import numpy as np
        points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        # preprocess.set_image("T1w", "/path/to/T1w.nii.gz")
        # preprocess.map_to_template(template_path = "/path/to/MNI152NLin2009bAsym.nii.gz", template_name = "MNI152NLin2009bAsym")
        mapped_points = preprocess.map_points_to_template(points, template_name = "MNI152NLin2009bAsym")
        ```

        '''
        map_info = self.get_template_mapping(template_name = template_name, native_type = native_type)
        if map_info is None:
            raise FileNotFoundError(f"Missing mapping from {native_type} to {template_name}. Please register the image to the template first.")
        map_args = map_info['native_to_template']
        import pandas as pd
        df = pd.DataFrame({
            'x' : - points[:, 0],
            'y' : - points[:, 1],
            'z' : points[:, 2]
        })
        mapped_points = ants.apply_transforms_to_points(
            dim = 3, points = df, 
            verbose = verbose,
            **map_args
        )
        return np.array([
            - mapped_points['x'].values,
            - mapped_points['y'].values,
            mapped_points['z'].values
        ]).transpose()
    
    def generate_atlas_from_template(self, template_name : str, template_atlas_folder : str, native_folder : Union[str, None] = None, verbose : bool=True) -> str:
        '''
        Generate a native atlas from the template atlas.

        @param template_name: The name of the template in which the atlases are defined, for example, "MNI152NLin2009bAsym"
        @type template_name: str

        @param template_atlas_folder: The path to the atlas folder or file (NIfTI).
        @type template_atlas_folder: str

        @param native_folder: The native folder name where the atlas will be saved. Default is inferred from `template_atlas_folder`.
        @type native_folder: str

        @param verbose: If True, print verbose output.
        @type verbose: bool

        @return: The path to the native atlas.
        @rtype: str
        '''
        # template_atlas_folder = "/Users/dipterix/Dropbox (PennNeurosurgery)/RAVE/Samples/Tower-related/templates 2/Lead-DBS_atlases_MNI_ICBM_2009b_NLIN_ASYM/MNI_ICBM_2009b_NLIN_ASYM/atlases/OCD Response Tract Atlas (Li 2020)"
        if not os.path.exists(template_atlas_folder) or not os.path.isdir(template_atlas_folder):
            raise FileNotFoundError(f"Invalid template atlas folder: {template_atlas_folder}")
        if native_folder is None:
            native_folder = os.path.basename(template_atlas_folder)
        for root_dir, dirnames, filenames in os.walk(template_atlas_folder):
            root_relpath = os.path.relpath(root_dir, template_atlas_folder)
            native_abspath = ensure_dir(file_path(self._work_path, "atlases", root_relpath))
            for filename in filenames:
                if filename.endswith(".nii.gz") and not filename.startswith("."):
                    template_path = file_path(root_dir, filename)
                    native_path = file_path(native_abspath, filename)
                    original_img = ants.image_read(template_path)
                    interpolator = "nearestNeighbor"
                    if original_img.max() <= 1.0 and original_img.min() >= 0.0:
                        if np.max(original_img[original_img < 1.0]) > 0.00001:
                            interpolator = "linear"
                    mapped_img = self.transform_image_from_template(
                        path = template_path, template_name = template_name,
                        native_type = "T1w", interpolator = interpolator, verbose = verbose)
                    mapped_img.to_file(native_path)
        return native_folder
        
