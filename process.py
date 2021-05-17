import os
import numpy as np
import nibabel as nib
from scipy.io import savemat, loadmat
import data_utils as du

try:
    import matlab.engine
    matlab_engine = matlab.engine.start_matlab()
    NO_ENGINE = False
except ImportError:
    print(
        "\n\n\nMatlab Engine API not found. Please execute the matlab script manually when prompted\n\n\n"
    )
    NO_ENGINE = True


def call_tapas_rDCM(header, time_series, region_labels, meta, manual):
    '''
    Parameters
    ----------
    header : nibabel.nifti1.Nifti1Header
        Nifti1 image header.
    time_series : numpy.ndarray
        2-dimensional array of shape (time_steps, regions) containg BOLD signal
        data.
    region_labels: List of strings
        Contains the names of regions.
    meta: dict
        Meta information to save alongside the data.

    Returns
    -------
    rDCK: some weird type #TODO: Figure this out.
        rDCM result.

    '''
    try:
        os.mkdir('.temp')
    except FileExistsError:
        pass

    if meta is None:
        meta = {'name': 'out'}

    # TODO: Make use of region names.
    to_mat = {
        'meta': meta,
        'Y': {
            'y': time_series,
            'dt': float(header.get_zooms()
                        [-1]),  # last value is time between scans in ms
            'name': region_labels
        }
    }

    savemat('.temp/in.mat', to_mat)

    # memory management to prevent me from dropping out of calls while running
    # the script in the future:
    del to_mat

    if manual or NO_ENGINE:
        print("Please execute the matlab script now.")
        print("After successful execution, press Enter to continue...")
        input()
    else:
        # call matlab script
        matlab_engine.call_tapas_rDCM(nargout=0)

    # not necessary since we're doing analysis in matlab:

    # rDCM = loadmat(f'.temp/{meta['name']}.mat')
    # print("\nSuccesfully loaded rDCM output.")
    # # TODO: figure out whether to delete .temp/
    # return rDCM


def rDCM_from_fMRI(filepaths, meta=None, scheme='yeo', manual=False):
    '''
    Parameters
    ----------
    filepaths : list
        List of paths pointing to .nii or .nii.gz files to use for rDCM.
    scheme: String or tuple
        Scheme to use for parcellation.

    Returns
    -------
    rDCK: some weird type #TODO: Figure this out.
        rDCM result.

    '''
    images = []

    for path in filepaths:
        images.append(nib.load(path))

    header0 = images[0].header

    if type(scheme) == str:  # if scheme isn't a masker already, make it one.
        scheme = du.make_masker(scheme)
    masker, region_labels = scheme

    # TODO: Make use of region names.
    time_series = du.parcellation(masker, images[0])
    for img in images[1:]:
        pi = du.parcellation(masker, img)
        time_series = du.combine_series(time_series, pi)

    return call_tapas_rDCM(header0, time_series, region_labels, meta, manual)


def batch_rDCM_from_fMRI(scheme_name, manual, data):
    if manual == True:
        print(
            "Warning: Using batch and manual mode simultaneously is antithetical."
        )

    scheme = du.make_masker(scheme_name)
    for key in data:
        ic(data[key])
        meta = {
            'name': key,
            'id': data[key]['id'],
            'scheme': scheme_name,
            'type': data[key]['type'],
        }
        rDCM_from_fMRI(data[key]['images'],
                       meta=meta,
                       scheme=scheme,
                       manual=manual)


if __name__ == '__main__':
    # example call
    rDCM_from_fMRI(
        ['path/to/file01.nii.gz', 'path/to/file02.nii.gz'],
        scheme='yeo',
    )
