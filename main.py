"""Expected function call:
python main.py path/to/file01.nii.gz [path/to/file02.nii.gz ...] [--scheme (harvox | yeo | aal)] [--manual]
"""
import os
import argparse
import numpy as np
import nibabel as nib
from scipy.io import savemat, loadmat
import data_utils as du

parser = argparse.ArgumentParser()
parser.add_argument('filepaths',
                    help='Path(s) to .nii or .nii.gz file(s).',
                    nargs='+')
parser.add_argument('--scheme',
                    default='yeo',
                    help='Which parcellation scheme to use.')
parser.add_argument('--manual',
                    dest='manual',
                    action='store_true',
                    help='Run matlab code manually.')
parser.set_defaults(manual=False)

args = parser.parse_args()

FILEPATHS = args.filepaths  #.split()
SCHEME = args.scheme
MANUAL_MATLAB = args.manual

try:
    import matlab.engine
    matlab_engine = matlab.engine.start_matlab()
except ImportError:
    print(
        "\n\n\nMatlab Engine API not found. Please execute the matlab script manually when prompted\n\n\n"
    )
    MANUAL_MATLAB = True


def call_tapas_rDCM(header, time_series):
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

    if MANUAL_MATLAB:
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


def rDCM_from_fMRI(filepaths):
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
    time_series = du.parcellation(SCHEME, images[0])[0]

    for img in images[1:]:
        pi = du.parcellation(SCHEME, img)[0]
        time_series = du.combine_series(time_series, pi)

    return call_tapas_rDCM(header0, time_series)


if __name__ == '__main__':
    rDCM_from_fMRI(FILEPATHS)
