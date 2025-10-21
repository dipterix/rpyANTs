#########################################################
# Custom Normalization scripts
#########################################################

import os
import ants
import numpy as np
from typing import Union
from ..utils.paths import normalize_path, file_path, try_import_antspynet
from ..utils.internals import get_lib_fn, ants_process_arguments
from ..utils.cache import _, StageContext, stage_image, restore_image

if False:
    import os
    import ants
    from typing import Union
    from rpyants.utils.paths import file_path, ensure_basename, file_copy, file_move, normalize_path, unlink, try_import_antspynet
    from rpyants.utils.internals import get_lib_fn, ants_process_arguments
    from rpyants.utils.cache import _, StageContext, stage_image, restore_image
    from rpyants.registration.normalization import normalization_with_atropos
    working_path = "/Users/dipterix/rave_data/raw_dir/Liming03/rave-imaging/work"
    mov_paths = ["/Users/dipterix/rave_data/raw_dir/Liming03/rave-imaging/3d.nii.gz"]
    fix_path = "/Users/dipterix/Library/Application Support/org.R-project.R/R/rpyANTs/templates/mni_icbm152_nlin_asym_09b/T1.nii.gz"
    weights = 1
    verbose = True
    use_antspynet = True
if False:
    normalization_with_atropos(
        fix_path = "/Users/dipterix/rave_data/others/three_brain/templates/mni_icbm152_nlin_asym_09b/T1.nii.gz",
        mov_paths = ["/Users/dipterix/rave_data/raw_dir/PAV038/rave-imaging/inputs/MRI/MRI_RAW.nii.gz"],
        working_path = "/Users/dipterix/rave_data/raw_dir/PAV038/rave-imaging/work",
        weights = 1,
        use_antspynet = True,
        verbose = True
    )

def normalize_to_template_syn(
    fix_path: Union[list, tuple, str], mov_paths: Union[list, tuple], outprefix: str, 
    weights: Union[float, int, list, tuple] = 1, verbose: bool = True,
    rigid_parameters = {
        "convergence": "[1000x500x250x0,1e-7,10]",
        "shrink-factors": "12x8x4x1",
        "smoothing-sigmas": "5x4x3x1vox",
    },
    affine_parameters = {
        "convergence": "[1000x500x250x0,1e-7,10]",
        "shrink-factors": "12x8x4x1",
        "smoothing-sigmas": "5x4x3x1vox",
    },
    syn_parameters = {
        "convergence": "[1000x500x500x0,1e-7,7]",
        "shrink-factors": "8x4x4x1",
        "smoothing-sigmas": "4x3x1x1vox",
    },
    refine_parameters = {
        "convergence": "[200x50x10x0,1e-6,7]",
        "shrink-factors": "4x4x2x1",
        "smoothing-sigmas": "2x2x1x1vox",
    }):
    '''
    Normalize images to template via volumetric mapping in 4 stages
    
    @param fix_path: Absolute path(s) to the fixing (template) image(s); can be a list or tuple, but must have the same lengths as `mov_paths` if do
    @type fix_path: list | tuple | str
    
    @param mov_paths: A list or tuple of absolute paths to the moving (native) images; all images should be registered to the same space
    @type mov_paths: list | tuple
    
    @param outprefix: Output path prefix
    @type outprefix: str
    
    @param weights: A single number or a list of numbers indicating the weights for each moving image
    @type weights: float | int | list | tuple
    
    @param rigid_parameters: dictionary of stage-1 rigid parameters
    @type rigid_parameters: dict; default is 
    {
        "convergence": "[1000x500x250x0,1e-7,10]",
        "shrink-factors": "12x8x4x1",
        "smoothing-sigmas": "5x4x3x1vox",
    }
    
    @param affine_parameters: dictionary of stage-2 affine parameters
    @type affine_parameters: dict; default is
    {
        "convergence": "[1000x500x250x0,1e-7,10]",
        "shrink-factors": "12x8x4x1",
        "smoothing-sigmas": "5x4x3x1vox",
    }
    
    @param syn_parameters: dictionary of stage-3 SyN parameters
    @type syn_parameters: dict; default is
    {
        "convergence": "[1000x500x500x0,1e-7,7]",
        "shrink-factors": "8x4x4x1",
        "smoothing-sigmas": "4x3x1x1vox",
    }
    
    @param refine_parameters: dictionary of stage-4 refine (also SyN) parameters
    @type refine_parameters: dict; default is
    {
        "convergence": "[200x50x10x0,1e-6,7]",
        "shrink-factors": "4x4x2x1",
        "smoothing-sigmas": "2x2x1x1vox",
    }
    
    @param verbose: Whether to verbose outputs
    @type verbose: bool
    '''
    # use forward slash to avoid escaping issues on windows
    # Do not use memory pointers!!!
    if isinstance(fix_path, str):
        fix_path = [fix_path]
    fix_path = [normalize_path(x, sep="/") for x in fix_path]
    if not isinstance(mov_paths, list):
        raise TypeError("normalize_to_template_syn: `mov_paths` must be a list")
    mov_paths = [normalize_path(moving_path, sep="/") for moving_path in mov_paths]
    if len(fix_path) != len(mov_paths):
        if len(fix_path) == 1:
            fix_path = fix_path * len(mov_paths)
        else:
            raise ValueError("normalize_to_template_syn: `fix_path` and `mov_paths` must have the same length or `fix_path` must be a single path")
    outprefix = normalize_path(outprefix, sep="/")
    verbose_flag = "1" if verbose else "0"
    if isinstance(weights, (float, int)):
        weights = [weights] * len(mov_paths)
    elif len(weights) == 1:
        weights = [weights[0]] * len(mov_paths)
    if len(weights) != len(mov_paths):
        raise ValueError("normalize_to_template_syn: weights must be a single number or a list of the same length as `mov_paths`")
    libfn = get_lib_fn("antsRegistration")
    # Build arguments, see https://github.com/ANTsX/ANTsPy/blob/78dad33ca4e12ae605e80d8d990ec73a52abe273/ants/registration/interface.py
    metrics = []
    for fix_path0, moving_path, weight in zip(fix_path, mov_paths, weights):
        metrics.append("--metric")
        metrics.append("MI[%s,%s,%.6f,32,Random,0.25]" % (fix_path0, moving_path, weight))
    args = [
        "--dimensionality", "3",
        "--float", "1",
        "-o", "[%s,%s,%s]" % (outprefix, outprefix + "Warped.nii.gz", outprefix + "InverseWarped.nii.gz"),
        "--interpolation", "Linear",
        "--use-histogram-matching", "1",
        "-v", verbose_flag,
        # Stage 1: Rigid Transformation
        "--transform", "Rigid[2.0]",
    ] + metrics + [
        "--convergence", rigid_parameters.get("convergence", "[1000x500x250x0,1e-7,10]"),
        "--shrink-factors", rigid_parameters.get("shrink-factors", "12x8x4x1"),
        "--smoothing-sigmas", rigid_parameters.get("smoothing-sigmas", "5x4x3x1vox"),
        # Stage 2: Affine Transformation
        "--transform", "Affine[0.1]",
    ] + metrics + [
        "--convergence", affine_parameters.get("convergence", "[1000x500x250x0,1e-7,10]"),
        "--shrink-factors", affine_parameters.get("shrink-factors", "12x8x4x1"),
        "--smoothing-sigmas", affine_parameters.get("smoothing-sigmas", "5x4x3x1vox"),
        # Stage 3: SyN Transformation
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
    if verbose_flag == "1":
        print("---- Call arguments --------------------------------")
        print("antsRegistration.sh %s" % " ".join(processed_args))
    if not reg_exit == 0:
        print(f"Registration failed with error code {reg_exit}")
    return {
        "outprefix": outprefix,
        "warpedmovout": outprefix + "Warped.nii.gz",
        "warpedfixout": outprefix + "InverseWarped.nii.gz",
        "fwdtransforms": [
            outprefix + "1Warp.nii.gz",
            outprefix + "0GenericAffine.mat", 
        ],
        "invtransforms": [
            outprefix + "0GenericAffine.mat", 
            outprefix + "1InverseWarp.nii.gz"
        ]
    }

def normalization_with_atropos(
    fix_path: Union[list, tuple, str], mov_paths: Union[list, tuple], working_path: str, 
    weights: Union[float, int, list, tuple] = 1, with_skull = False, cost_function = "CC",
    use_antspynet: bool = True, 
    verbose: bool = True) -> dict:
    '''
    Normalize images to template with ANTs Atropos as constraints.

    @param fix_path: Absolute path(s) to the fixing (template) image(s); can be a list or tuple, but must have the same lengths as `mov_paths` if do
    @type fix_path: list | tuple | str

    @param mov_paths: A list or tuple of absolute paths to the moving (native) images; all images should be registered to the same space
    @type mov_paths: list | tuple

    @param working_path: Working directory path
    @type working_path: str

    @param weights: A single number or a list of numbers indicating the weights for each moving image
    @type weights: float | int | list | tuple

    @param with_skull: Whether to normalize with skulls or skull-striped
    @type with_skull: bool

    @param cost_function: Cost function for the T1w image; default is 'CC'; choices are 'MI', 'CC'
    @type cost_function: str

    @param use_antspynet: Whether to use antspynet for deep_atropos
    @type use_antspynet: bool

    @param verbose: Whether to verbose outputs; default is True
    @type verbose: bool

    @return: Registration results
    '''
    # Check whether antspynet is available
    if use_antspynet:
        antspynet = try_import_antspynet()
    else:
        antspynet = None
    # make fix_path a list
    if isinstance(fix_path, str):
        fix_paths = [fix_path]
    else:
        fix_paths = fix_path
    fix_paths = [normalize_path(x, sep="/") for x in fix_paths]
    # make mov_paths a list
    if not isinstance(mov_paths, (list, tuple)):
        mov_paths = [mov_paths]
    mov_paths = [normalize_path(moving_path, sep="/") for moving_path in mov_paths]

    n_images = len(mov_paths)
    if n_images == 0:
        raise ValueError("normalization_with_atropos: `mov_paths` must have at least one image")
    if n_images > 1:
        if len(fix_paths) == 1:
            fix_paths = fix_paths * n_images
        elif len(fix_paths) != n_images:
            raise ValueError("normalization_with_atropos: `fix_path` and `mov_paths` must have the same length or `fix_path` must be a single path")
    if isinstance(weights, (float, int)):
        weights = [weights] * n_images
    if len(weights) != n_images:
        if len(weights) == 1:
            weights = [weights[0]] * n_images
        else:
            raise ValueError("normalization_with_atropos: weights must be a single number or a list of the same length as `mov_paths`")
    
    # Treats the first image as the T1 image
    fixing_path = fix_paths[0]
    fixing_root = os.path.dirname(fixing_path)
    fixing_mask_path = file_path(fixing_root, "T1_brainmask.nii.gz")
    moving_path = mov_paths[0]
    if not isinstance(weights, list):
        weights = [weights for x in mov_paths]
    if len(weights) != len(mov_paths):
        raise ValueError("normalization_with_atropos: weights must be a single number or a list of the same length as `mov_paths`")
    # Read in image 
    with StageContext("T1", "image", working_path, verbose=verbose) as (ctx, moving_img):
        if moving_img is None:
            moving_img = ants.image_read(moving_path)
            ctx.result = moving_img
    
    # apply intensity normalization (pass 1)
    with StageContext("nu", "image", working_path, verbose=verbose) as (ctx, _._):
        if ctx.result is None:
            ctx.result = ants.abp_n4(moving_img)
        nu_img = ctx.result
    # # apply intensity normalization
    # with StageContext("nu3", "image", working_path, verbose=verbose) as (ctx, _._):
    #     if ctx.result is None:
    #         ctx.result = ants.abp_n4(moving_img, usen3=True)
    fixing_img = ants.image_read(fixing_path)
    fixing_mask = ants.image_read(fixing_mask_path)

    stage_image(fixing_img, name="fixing.nii.gz", root=working_path)

    # n4_bias_field_correction(image, mask=None, rescale_intensities=False, shrink_factor=4, convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07}, spline_param=None, return_bias_field=False, verbose=False, weight_mask=None)
    # antspynet is not available, use ants.atropos
    # Need to get brain mask and find initial registration to template
    # Register brain using `SyNabp`: SyN optimized for abpBrainExtraction.
    with StageContext("SyNabp", "registration", working_path) as (ctx, res_abp):
        if res_abp is None:
            res_abp = ants.registration(
                fixed=fixing_img, moving=nu_img,
                type_of_transform='SyNabp', verbose=verbose)
            ctx.result = res_abp
    
    # Quick 3-way atropos to get white-matter, or whatever the brightest
    # Use the mask to normalize intensity again
    with StageContext("nu2", "image", working_path, verbose=verbose) as (ctx, _._):
        if ctx.result is None:
            # Quick brain-mask
            brain_mask0 = ants.apply_transforms(
                fixed=moving_img, 
                moving=fixing_mask, 
                transformlist=res_abp["invtransforms"], 
                interpolator="nearestNeighbor")
            # Quick 3-way atropos to get white-matter, or whatever the brightest
            seg3 = ants.atropos(a=moving_img, x = brain_mask0, i="kmeans[3]", m="[0.1,1x1x1]", c="[5,0]")
            probs3 = seg3["probabilityimages"]  # list of 3 probability maps
            means = []
            for p in probs3:
                mu = (moving_img*p).sum() / (p.sum() + 1e-8)
                means.append(mu)
            wm_idx = int(np.argmax(means))
            wm_mean = float(means[wm_idx])
            scale = 110.0 / (wm_mean if wm_mean > 0 else 1.0)
            t1_norm = moving_img * scale                   # like FS T1.mgz inside brain
            ctx.result = ants.abp_n4(moving_img)
            wm_prob = probs3[wm_idx]
            n4c = ants.n4_bias_field_correction(
                t1_norm, weight_mask=wm_prob,
                rescale_intensities=True,
                shrink_factor=1,
                convergence={'iters':[50,50,30,20],'tol':1e-7}
                # spline_param=200
            )
            ctx.result = n4c
        moving_img = ctx.result

    # Denoise image and also smooth intensities
    # dn = ants.denoise_image(moving_img, shrink_factor=1, p=1, r=3, noise_model="Rician")

    # Atropos: Check whether antspynet is available
    if antspynet is None:
        
        # extract brain 
        with StageContext("brainmask", "image", working_path) as (ctx, brain_mask):
            if brain_mask is None:
                brain_mask = ants.apply_transforms(
                    fixed=moving_img, 
                    moving=fixing_mask, 
                    transformlist=res_abp["invtransforms"], 
                    interpolator="nearestNeighbor")
                stage_image(moving_img * brain_mask, "brain.nii.gz", root=working_path)
                ctx.result = brain_mask

        with StageContext("atropos", "image[5]", working_path) as (ctx, moving_atropos):
            if moving_atropos is None:
                # Atropos to 6+1 segments (CSF + GM + WM + ...)
                with StageContext("atropos_prior", "image[6]", working_path) as (ctx_prior, moving_atropos_prior):
                    if moving_atropos_prior is None:
                        moving_atropos_prior = []
                        for ii in range(6):
                            moving_atropos_prior.append(
                                ants.apply_transforms(
                                    fixed=moving_img, 
                                    moving=ants.image_read(file_path(fixing_root, f"atropos_{ii + 2}.nii.gz")), 
                                    transformlist=res_abp["invtransforms"], 
                                    interpolator="linear",
                                    verbose=verbose
                                )
                            )
                        ctx_prior.result = moving_atropos_prior
                ctx_prior = None
                moving_atropos_prior2 = [
                    moving_atropos_prior[0] * 0.2 + 0.4, # CSF from 0.4 - 0.6 prob
                    moving_atropos_prior[1] * 0.6 + 0.2, # GM: from 0.2 - 0.8
                    (moving_atropos_prior[2] + moving_atropos_prior[3] + moving_atropos_prior[4]) * 0.6 + 0.2, # WM, Thalamus, Brain Stem
                    moving_atropos_prior[5], # Cerebellum
                ]
                moving_atropos_prior = None
                moving_atropos_res = ants.atropos(
                    a=moving_img, m='[0.1,1x1x1]', c='[25,0]',
                    i=moving_atropos_prior2, x=brain_mask, priorweight=0.0, 
                    verbose=1 if verbose else 0)
                moving_atropos = moving_atropos_res['probabilityimages']
                moving_atropos.insert(0, moving_atropos_res['segmentation'])
                # StageContext("atropos", "image[5]", working_path)._store_result(moving_atropos)
                ctx.result = moving_atropos
                moving_atropos_res = None

        with StageContext("SyN_w_atropos", "registration", working_path) as (ctx, res_syn):
            if res_syn is None:
                fixing_csf = ants.image_read(file_path(fixing_root, "atropos_2.nii.gz"))
                # fixing_csf_path = file_path(fixing_root, "atropos_2.nii.gz")
                with StageContext("atropos_dGWbS", "image", fixing_root) as (ctx, fixing_dGWbS):
                    if fixing_dGWbS is None:
                        fixing_dGWbS = ants.image_read(file_path(fixing_root, "atropos_4.nii.gz")) + ants.image_read(file_path(fixing_root, "atropos_5.nii.gz")) + ants.image_read(file_path(fixing_root, "atropos_6.nii.gz"))
                        ctx.result = fixing_dGWbS
                # fixing_dGWbS_path = file_path(fixing_root, "atropos_dGWbS.nii.gz")
                moving_csf = moving_atropos[1]
                # moving_csf_path = file_path(working_path, "atropos_1.nii.gz")
                moving_dGWbS = moving_atropos[3]
                # moving_dGWbS_path = file_path(working_path, "atropos_3.nii.gz")
                # fix_paths2 = fix_paths + [fixing_csf_path, fixing_dGWbS_path]
                fix_paths2 = [x for x in fix_paths]
                fix_paths2.pop(0)
                fix_paths2 = [ants.image_read(x) for x in fix_paths2] + [fixing_csf, fixing_dGWbS]

                # mov_paths2 = mov_paths + [moving_csf_path, moving_dGWbS_path]
                mov_paths2 = [x for x in mov_paths]
                mov_paths2.pop(0)
                mov_paths2 = [ants.image_read(x) for x in mov_paths2] + [moving_csf, moving_dGWbS]
                
                weights2 = weights + [0.1, 0.5]
                weights2.pop(0)
                multivariate_extras = [
                    ("MI", fpath, mpath, weight, '32,Random,0.25') for fpath, mpath, weight in zip(fix_paths2, mov_paths2, weights2)
                ]
                stage_image(fixing_img * fixing_mask, "fixed_brain.nii.gz", root=working_path)
                mov_skullstrip = restore_image("brain.nii.gz", root=working_path)
                fix_skullstrip = restore_image("fixed_brain.nii.gz", root=working_path)
                # res_syn = ants.registration(
                #     fixed=fixing_img,
                #     moving=moving_img,
                #     type_of_transform='antsRegistrationSyNQuick[so]',
                #     mask=fixing_mask,
                #     moving_mask=brain_mask,
                #     mask_all_stages=False,
                #     initial_transform=[res_abp['fwdtransforms'][1]],
                #     multivariate_extras=multivariate_extras, 
                #     verbose=verbose
                # )
                syn_sampling = 32
                if cost_function == 'CC':
                    syn_sampling = 4   # CC radius ~4
                if with_skull:
                    res_syn = ants.registration(
                        fixed=fixing_img,
                        moving=moving_img,
                        type_of_transform='SyNAggro',
                        flow_sigma=3.5, total_sigma=0,        # mild extra smoothing of total field
                        aff_metric='mattes', aff_sampling=32, 
                        aff_random_sampling_rate=0.2, 
                        syn_metric=cost_function, syn_sampling=syn_sampling,
                        reg_iterations=(100,70,50,0),    
                        mask=fixing_mask,
                        moving_mask=brain_mask,
                        mask_all_stages=True,
                        multivariate_extras=multivariate_extras, 
                        verbose=verbose
                    )
                else:
                    res_syn = ants.registration(
                        fixed=fix_skullstrip,
                        moving=mov_skullstrip,
                        type_of_transform='SyNAggro',
                        grad_step = 0.15, 
                        flow_sigma=3.5, total_sigma=0,        # mild extra smoothing of total field
                        aff_metric='mattes', aff_sampling=32, 
                        aff_random_sampling_rate=0.2, 
                        syn_metric=cost_function, syn_sampling=syn_sampling,
                        reg_iterations=(100,70,50,0),    
                        mask=fixing_mask,
                        moving_mask=brain_mask,
                        mask_all_stages=True,
                        multivariate_extras=multivariate_extras, 
                        verbose=verbose
                    )
                # StageContext("SyN_w_atropos", "registration", working_path)._store_result(res_syn)
                ctx.result = res_syn
    else:
        # antspynet is available, use deep_atropos
        # with StageContext("atropos", "image[8]", fixing_root) as (ctx, fixing_deep_atropos_list):
        #     if fixing_deep_atropos_list is None:
        #         fixing_deep_atropos = antspynet.deep_atropos(
        #             t1=fixing_img,
        #             do_preprocessing=True,
        #             use_spatial_priors=1,
        #             verbose=verbose)
        #         fixing_deep_atropos_list = fixing_deep_atropos['probability_images']
        #         fixing_deep_atropos_list.insert(0, fixing_deep_atropos['segmentation_image'])
        #         fixing_deep_atropos = None
        #         ctx.result = fixing_deep_atropos_list

        with StageContext("deep_atropos", "image[8]", working_path) as (ctx, moving_deep_atropos_list):
            if moving_deep_atropos_list is None:
                moving_deep_atropos = antspynet.deep_atropos(
                    t1=moving_img,
                    do_preprocessing=True,
                    use_spatial_priors=1,
                    verbose=verbose)
                moving_deep_atropos_list = moving_deep_atropos['probability_images']
                moving_deep_atropos_list.insert(0, moving_deep_atropos['segmentation_image'])
                moving_deep_atropos = None
                ctx.result = moving_deep_atropos_list
        
        # extract brain 
        with StageContext("brainmask", "image", working_path) as (ctx, brain_mask):
            if brain_mask is None:
                brain_mask = moving_deep_atropos_list[0] > 0
                # brain_mask = ants.apply_transforms(
                #     fixed=moving_img, 
                #     moving=fixing_mask, 
                #     transformlist=res_abp["invtransforms"], 
                #     interpolator="nearestNeighbor")
                stage_image(moving_img * brain_mask, "brain.nii.gz", root=working_path)
                ctx.result = brain_mask
        
        # extracted brain
        with StageContext("brain", "image", working_path) as (ctx, brain_skulstrip):
            if brain_skulstrip is None:
                brain_skulstrip = moving_img * brain_mask
                ctx.result = brain_skulstrip
        
        # fixing_deep_atropos_list = None
        stage_image(fixing_img * fixing_mask, "fixed_brain.nii.gz", root=working_path)
        brain_mask_no_csf = moving_deep_atropos_list[0] > 0
        brain_mask_no_csf[moving_deep_atropos_list[1] > 0.05] = 0
        stage_image(moving_img * brain_mask_no_csf, "brain.nii.gz", root=working_path)
        mov_skullstrip = restore_image("brain.nii.gz", root=working_path)
        fix_skullstrip = restore_image("fixed_brain.nii.gz", root=working_path)
        moving_deep_atropos_list = None
        brain_mask_no_csf = None

        # gc.collect()
        with StageContext("SyN_w_deep_atropos", "registration", working_path) as (ctx, res_syn):
            if res_syn is None:
                # fixing_csf = fixing_deep_atropos_list[2]
                fixing_csf_path = file_path(fixing_root, "atropos_2.nii.gz")
                # fixing_deepGM = fixing_deep_atropos_list[5]
                fixing_deepGM_path = file_path(fixing_root, "atropos_5.nii.gz")
                # moving_csf = moving_deep_atropos_list[2]
                moving_csf_path = file_path(working_path, "deep_atropos_2.nii.gz")
                # moving_deepGM = moving_deep_atropos_list[5]
                moving_deepGM_path = file_path(working_path, "deep_atropos_5.nii.gz")
                fix_paths2 = fix_paths + [fixing_csf_path, fixing_deepGM_path]
                mov_paths2 = mov_paths + [moving_csf_path, moving_deepGM_path]
                weights2 = weights + [0.5, 0.5]
                fix_paths2.pop(0)
                mov_paths2.pop(0)
                weights2.pop(0)
                multivariate_extras = [
                    ("MI", ants.image_read(fpath), ants.image_read(mpath), weight, '32,Random,0.25') for fpath, mpath, weight in zip(fix_paths2, mov_paths2, weights2)
                ]
                # with StageContext("nu2", "image", working_path, verbose=verbose) as (ctx, _._):
                # registration(fixed, moving, type_of_transform='SyN', initial_transform=None, outprefix='', mask=None, moving_mask=None, mask_all_stages=False, 
                #              grad_step=0.2, flow_sigma=3, total_sigma=0, aff_metric='mattes', aff_sampling=32, aff_random_sampling_rate=0.2, 
                #              syn_metric='mattes', syn_sampling=32, reg_iterations=(40, 20, 0), aff_iterations=(2100, 1200, 1200, 10), aff_shrink_factors=(6, 4, 2, 1), 
                #              aff_smoothing_sigmas=(3, 2, 1, 0), write_composite_transform=False, random_seed=None, verbose=False, multivariate_extras=None, 
                #              restrict_transformation=None, smoothing_in_mm=False, singleprecision=True, use_legacy_histogram_matching=False, **kwargs)
                syn_sampling = 32
                if cost_function == 'CC':
                    syn_sampling = 4   # CC radius ~4
                if with_skull:
                    res_syn = ants.registration(
                        fixed=fixing_img,
                        moving=moving_img,
                        type_of_transform='SyNAggro',
                        flow_sigma=3.5, total_sigma=0,        # mild extra smoothing of total field
                        aff_metric='mattes', aff_sampling=32, 
                        aff_random_sampling_rate=0.2, 
                        syn_metric=cost_function, syn_sampling = syn_sampling,
                        reg_iterations=(100,70,50,0),    
                        mask=fixing_mask,
                        moving_mask=brain_mask,
                        mask_all_stages=True,
                        multivariate_extras=multivariate_extras, 
                        verbose=verbose
                    )
                else:
                    res_syn = ants.registration(
                        fixed=fix_skullstrip,
                        moving=mov_skullstrip,
                        type_of_transform='SyNAggro',
                        grad_step = 0.15, 
                        flow_sigma=3.5, total_sigma=0,        # mild extra smoothing of total field
                        aff_metric='mattes', aff_sampling=32, 
                        aff_random_sampling_rate=0.2, 
                        syn_metric=cost_function, syn_sampling=syn_sampling,        
                        reg_iterations=(100,70,50,0),    
                        mask=fixing_mask,
                        moving_mask=brain_mask,
                        mask_all_stages=True,
                        multivariate_extras=multivariate_extras, 
                        verbose=verbose
                    )
                # initial_transform=[res_abp['fwdtransforms'][1]],
                ctx.result = res_syn
                # StageContext("SyN_w_deep_atropos", "registration", working_path)._store_result(res_syn)
    return res_syn
