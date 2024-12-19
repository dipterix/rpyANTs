import numpy as np

M44_LPS_to_RAS = np.array([
    [-1., 0., 0., 0.],
    [0., -1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 1.],
])

def ants_AffineTransform_to_m44( affine ):
    x_new = affine.apply_to_vector([1., 0., 0.])
    y_new = affine.apply_to_vector([0., 1., 0.])
    z_new = affine.apply_to_vector([0., 0., 1.])
    origin = affine.apply_to_point([0., 0., 0.])
    m34 = np.array([x_new, y_new, z_new, origin]).transpose()
    # Append the last row
    m44 = np.vstack([m34, [0., 0., 0., 1.]])
    return m44

def get_xform(img):
    '''
    Get `xform` (Voxel to RAS) transform from ANTsImage
    @param img: an ANTsImage instance
    @type img: ANTsImage
    '''
    ndims = img.dimension
    # each column of the direction is multiplied by an elem of spacing
    # to rescale
    # np.matmul(img.direction, np.diag(spacing))
    linear_matrix = img.direction * img.spacing
    xform_lps = np.vstack(
        (
            np.hstack((linear_matrix, np.array(img.origin).reshape(ndims, 1))),
            np.array([1. if x == ndims else 0. for x in range(ndims + 1)]).reshape(1, ndims + 1)
        )
    )
    lps_to_ras = np.diag([-1. if x < 2 else 1. for x in range(ndims + 1)])
    return np.matmul(lps_to_ras, xform_lps)

def set_xform(img, vox2ras):
    # import rpyants
    # vox2ras = rpyants.utils.transforms.get_xform(img)
    ndims = img.dimension
    ras_to_lps = np.diag([-1. if x < 2 else 1. for x in range(ndims + 1)])
    vox2lps = np.matmul(ras_to_lps, vox2ras)
    # get translate
    translate = vox2lps[:ndims, ndims]
    # linear part
    linear_matrix = vox2lps[:ndims, :ndims]
    # normalize each column, axis=0 means L2 norm for each column
    spacing = np.sqrt(np.sum(np.power(linear_matrix, 2.), 0))
    direction = linear_matrix / spacing
    img.set_direction(direction)
    img.set_spacing(spacing.tolist())
    img.set_origin(translate.tolist())
    return img

def apply_transform_to_nifti(src_path, dst_path, ras_to_ras):
    # nii = nib.load(src_path)
    # new_sform = np.matmul(ras_to_ras, nii.get_sform())
    # nii.set_sform(new_sform)
    # nii.set_qform(new_sform)
    # nib.save(nii, dst_path)
    import ants
    nii = ants.image_read(src_path)
    new_sform = np.matmul(ras_to_ras, get_xform(nii))
    nii = set_xform(nii, new_sform)
    nii.to_file(dst_path)
    return dst_path
