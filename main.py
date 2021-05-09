import os
import numpy as np
import nibabel as nib
from scipy.io import savemat, loadmat
import TNM_project_input_data as pid
import matlab.engine

matlab_engine = matlab.engine.start_matlab()


def call_tapas_rDCM(header, freq_series):
    try:
        os.mkdir('.temp')
    except FileExistsError:
        pass

    to_mat = {
        'a': [],
        'b': [],
        'c': [],
        'd': [],
        'Y': {
            'y': freq_series,
            'dt': 0.5
            # header.get_zooms()[-1]  # last value is time between scans in ms
        },
        'U': {
            'u': np.zeros((16 * freq_series.shape[0], 1)),
            # 'dt': header.get_zooms()[-1] / 16
        },
        'v': freq_series.shape[0],  # specify number of datapoints (per region)
        'n': freq_series.shape[1]  # specify number of regions
    }

    print(to_mat['Y']['y'].shape)
    print(to_mat['U']['u'].shape)

    savemat('.temp/in.mat', to_mat)

    # call matlab script
    matlab_engine.call_tapas_rDCM(nargout=0)

    rDCM = loadmat('.temp/out.mat')
    # add cleanup such that .temp/ is deleted?
    return rDCM


def rDCM_from_fMRI(path):
    img = nib.load(path)
    parcellated_img = pid.parcellation("Yeo", img)
    freq_series = pid.Fourier_transform(parcellated_img[0])
    return call_tapas_rDCM(img.header, freq_series)


def example_call():
    rDCM_from_fMRI(
        'data/ds003059/sub-001/ses-LSD/func/sub-001_ses-LSD_task-rest_run-01_bold.nii.gz'
    )


if __name__ == '__main__':
    example_call()