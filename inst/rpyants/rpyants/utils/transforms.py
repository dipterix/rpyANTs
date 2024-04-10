import numpy as np

def ants_AffineTransform_to_m44( affine ):
    x_new = affine.apply_to_vector([1., 0., 0.])
    y_new = affine.apply_to_vector([0., 1., 0.])
    z_new = affine.apply_to_vector([0., 0., 1.])
    origin = affine.apply_to_point([0., 0., 0.])
    m34 = np.array([x_new, y_new, z_new, origin]).transpose()
    # Append the last row
    m44 = np.vstack([m34, [0., 0., 0., 1.]])
    return m44

M44_LPS_to_RAS = np.array([
    [-1., 0., 0., 0.],
    [0., -1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 1.],
])

def apply_transform_to_nifti(src_path, dst_path, ras_to_ras):
    import nibabel as nib
    nii = nib.load(src_path)
    new_sform = np.matmul(ras_to_ras, nii.get_sform())
    nii.set_sform(new_sform)
    nii.set_qform(new_sform)
    nib.save(nii, dst_path)
    return dst_path

