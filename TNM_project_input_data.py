#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 17:27:57 2021

@author: zhenrujia
"""

import nibabel as nib
def load_data(path_name):
    '''
    Parameters
    ----------
    path_name : String
        Path of the data.

    Returns
    -------
    bold_data : Nifti1Image
    bold_data_df : memmap
        Extracted from Nifti1Image.
    '''
    bold_data = nib.load(path_name)
    bold_data_df = bold_data
    bold_data_df.get_fdata()
    bold_data_df = bold_data_df.get_fdata()
    print(bold_data_df.shape)
    return (bold_data,bold_data_df)




from nilearn import datasets
from nilearn.input_data import NiftiLabelsMasker

def parcellation(Type, data):
    '''
    Parameters
    ----------
    Type : String
        The type of parcellation wanted.
    data : Nifti1Image
        The data directly loaded, not the extracted one.

    Returns
    -------
    time_series: array
        Time series of the signal in each region (reduce the region with mean).
    labels: list
        Labels of all the regions in parcellation.
    '''
    if Type == "Harvard-Oxford": # 48 regions
        dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
        atlas_filename = dataset.maps
        labels = dataset.labels
        masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True, 
                                   high_variance_confounds = True, verbose = 1)
        time_series = masker.fit_transform(data)
        print ("Using the Harvard-Oxford parcellations, there are {} regions.".format(time_series.shape[1]))
    elif Type == "Yeo": # 17 regions
        dataset = datasets.fetch_atlas_yeo_2011()
        masker = NiftiLabelsMasker(labels_img=dataset['thick_17'], standardize=True, 
                                   high_variance_confounds = True, memory='nilearn_cache',  verbose = 1)
        time_series = masker.fit_transform(data)
        labels = list(np.arange(0,time_series.shape[1]+1))
        print ("Using the Yeo 2011 parcellation, there are {} regions.".format(time_series.shape[1]))

    elif Type == "AAL": # 116 regions
        dataset = datasets.fetch_atlas_aal(version='SPM12')
        labels = ["Background"] + dataset['labels']
        masker = NiftiLabelsMasker(labels_img=dataset['maps'], standardize=True, 
                                   high_variance_confounds = True, memory='nilearn_cache', verbose = 1)
        time_series = masker.fit_transform(data)
        print ("Using the AAL parcellation, there are {} regions.".format(time_series.shape[1]))
    print (time_series.shape)
    return time_series,labels




import numpy as np
def Fourier_transform(data):
    '''
    Parameters
    ----------
    data : array
        Time series data.

    Returns
    -------
    freq_series : array
        Discrete Fourier Transform into the frequency domain.
    '''
    freq_series = np.fft.fft(data, axis=0) #axis = 0: the transform is performed on each column
    return freq_series



from nilearn.connectome import ConnectivityMeasure
from nilearn import plotting
def plot_corr(time_series,labels):
    '''
    Parameters
    ----------
    time_series : array
        Time series of the signal in each region .
    labels : list
        Labels of all the regions in parcellation.

    Returns
    -------
    None.
    '''
    correlation_measure = ConnectivityMeasure(kind='correlation')
    correlation_matrix = correlation_measure.fit_transform([time_series])[0]

    # Mask the main diagonal for visualization:
    np.fill_diagonal(correlation_matrix, 0)
    # The labels we have start with the background (0), hence we skip the first label
    plotting.plot_matrix(correlation_matrix, figure=(10, 8), labels=labels[1:],
                         vmax=0.8, vmin=-0.8, reorder=True)





# =============================================================================
#  --- try to find Confound ---
# from nilearn.input_data import NiftiSpheresMasker
# seed_masker = NiftiSpheresMasker([(0, -53, 26)])
# seed_time_series = seed_masker.fit_transform(img)
# 
# t_r = 2.
# n_scans = 176
# frametimes = np.linspace(0, (n_scans - 1) * t_r, n_scans)
# 
# design_matrix = make_first_level_design_matrix(frametimes, hrf_model='spm',
#                                                add_regs=seed_time_series,
#                                                add_reg_names=["pcc_seed"])
# 
# first_level_model = FirstLevelModel(t_r=t_r, slice_time_ref=slice_time_ref)
# first_level_model = first_level_model.fit(run_imgs=adhd_dataset.func[0],
#                                     design_matrices=design_matrix)
# 
# =============================================================================
