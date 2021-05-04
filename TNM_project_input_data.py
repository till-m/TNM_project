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
        The type of parcellation wanted
    data : Nifti1Image
        The data directly loaded, not the extracted one.

    Returns
    -------
    time_series: array
        Time series of the signal in each region (need to read into how the signals are computed)
    '''
    if Type == "Harvard-Oxford": # 48 regions
        dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
        atlas_filename = dataset.maps
        labels = dataset.labels
        masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True)
        time_series = masker.fit_transform(data)
        print ("Using the Harvard-Oxford parcellations, there are {} regions.".format(time_series.shape[1]))
    elif Type == "Yeo": # 17 regions
        dataset = datasets.fetch_atlas_yeo_2011()
        masker = NiftiLabelsMasker(labels_img=dataset['thick_17'], standardize=True,
                           memory='nilearn_cache')
        time_series = masker.fit_transform(data)
        print ("Using the Yeo 2011 parcellation, there are {} regions.".format(time_series.shape[1]))

    elif Type == "AAL": # 116 regions
        dataset = datasets.fetch_atlas_aal(version='SPM12')
        masker = NiftiLabelsMasker(labels_img=dataset['maps'], standardize=True,
                           memory='nilearn_cache')
        time_series = masker.fit_transform(data)
        print ("Using the AAL parcellation, there are {} regions.".format(time_series.shape[1]))
    print (time_series.shape)
    return time_series




