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


def call_tapas_rDCM(header, time_series, manual):
    '''
    Parameters
    ----------
    header : nibabel.nifti1.Nifti1Header
        Nifti1 image header.
    time_series : numpy.ndarray
        2-dimensional array of shape (time_steps, regions) containg BOLD signal
        data.

    Returns
    -------
    rDCK: some weird type #TODO: Figure this out.
        rDCM result.

    '''
    try:
        os.mkdir('.temp')
    except FileExistsError:
        pass

    # TODO: Make use of region names.
    to_mat = {
        'y': time_series,
        'dt': 0.5
        # header.get_zooms()[-1]  # last value is time between scans in ms
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

    rDCM = loadmat('.temp/out.mat')
    print("\nSuccesfully loaded rDCM output.")
    # TODO: figure out whether to delete .temp/
    return rDCM


def rDCM_from_fMRI(filepaths, scheme='yeo', manual=False):
    '''
    Parameters
    ----------
    filepaths : list
        List of paths pointing to .nii or .nii.gz files to use for rDCM.

    Returns
    -------
    rDCK: some weird type #TODO: Figure this out.
        rDCM result.

    '''
    images = []
    for path in filepaths:
        images.append(nib.load(path))

    header0 = images[0].header

    # TODO: Make use of region names.
    time_series = du.parcellation(scheme, images[0])[0]

    for img in images[1:]:
        pi = du.parcellation(scheme, img)[0]
        time_series = du.combine_series(time_series, pi)

    return call_tapas_rDCM(header0, time_series, manual)


if __name__ == '__main__':
    # example call
    rDCM_from_fMRI(
        ['path/to/file01.nii.gz', 'path/to/file02.nii.gz'],
        scheme='yeo',
    )
