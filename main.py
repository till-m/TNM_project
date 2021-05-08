import os
import nibabel as nib
from scipy.io import savemat, loadmat
import TNM_project_input_data as pid
import matlab.engine

matlab_engine = matlab.engine.start_matlab()


def call_tapas_rDCM(freq_series):
    try:
        os.mkdir('./.temp')
    except FileExistsError:
        pass

    to_mat = {'Y': {'y': freq_series}}  # TODO: add the other fields

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
    return call_tapas_rDCM(freq_series)


def example_call():
    rDCM_from_fMRI(
        'data/ds003059/sub-001/ses-LSD/func/sub-001_ses-LSD_task-rest_run-01_bold.nii.gz'
    )


if __name__ == '__main__':
    example_call()