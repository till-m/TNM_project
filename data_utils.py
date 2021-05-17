#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 17:27:57 2021

@author: zhenrujia
"""

import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import datasets
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from nilearn import plotting


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
    bold_data_df = bold_data.get_fdata()
    print("The shape of the data is", bold_data_df.shape)
    return (bold_data, bold_data_df)


def make_masker(scheme):
    '''
    Parameters
    ----------
    scheme : String
        The type of parcellation wanted.

    Returns
    -------
    masker: nilearn.input_data.NiftiLabelsMasker
        Masker of the chosen scheme.
    labels: list
        Labels of all the regions in parcellation.
    '''
    if scheme.lower() == "harvox":  # 48 regions
        dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
        atlas_filename = dataset.maps
        labels = dataset.labels
        masker = NiftiLabelsMasker(labels_img=atlas_filename,
                                   standardize=True,
                                   high_variance_confounds=True,
                                   verbose=1)
    elif scheme.lower() == "yeo":  # 17 regions
        dataset = datasets.fetch_atlas_yeo_2011()
        masker = NiftiLabelsMasker(labels_img=dataset['thick_17'],
                                   standardize=True,
                                   high_variance_confounds=True,
                                   verbose=1)
        labels = list(np.arange(1, 18))
    elif scheme.lower() == "aal":  # 116 regions
        dataset = datasets.fetch_atlas_aal(version='SPM12')
        labels = ["Background"] + dataset['labels']
        masker = NiftiLabelsMasker(labels_img=dataset['maps'],
                                   standardize=True,
                                   high_variance_confounds=True,
                                   verbose=1)
    return masker, labels


def parcellation(masker, data):
    '''
    Parameters
    ----------
    masker: nilearn.input_data.NiftiLabelsMasker
        Masker to use.
    data : Nifti1Image
        The data directly loaded, not the extracted one.

    Returns
    -------
    time_series: array
        Time series of the signal in each region (reduce the region with mean).
    '''

    time_series = masker.fit_transform(data)
    print(
        f"Using the chosen parcellation scheme, there are {time_series.shape[1]} regions with {time_series.shape[0]} datapoints each."
    )
    return time_series


def combine_series(time_series_1, time_series_2):
    '''
    Parameters
    ----------
    time_series_1 : array
        Array of time series 01.
    time_series_2 : array
        Array of time series 03.

    Returns
    -------
    combined_series : array
        Combined array of time series 01 and time series 03.
    '''
    print("Combining the two time series")
    combined_series = np.concatenate((time_series_1, time_series_2), axis=0)
    return combined_series


def plot_corr(time_series, labels):
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
    plotting.plot_matrix(correlation_matrix,
                         figure=(10, 8),
                         labels=labels[1:],
                         vmax=0.8,
                         vmin=-0.8,
                         reorder=True)


def save_data(data, subject, LSD):
    '''
    Parameters
    ----------
    data : array
        Frequency series.
    subject : string
        Number of the subject.
    LSD: int
        Indicator of whether it's the data after taking LSD.

    Returns
    -------
    pathname : string
        Path name of the file saved.
    '''
    print("Saving the data")
    data_df = pd.DataFrame(data)
    if LSD == 1:
        pathname = '/Users/zhenrujia/Downloads/Time_series_BOLD' + '_sub' + subject + '_LSD' + '.csv'
    else:
        pathname = '/Users/zhenrujia/Downloads/Time_series_BOLD' + '_sub' + subject + '_PLCB' + '.csv'
    data_df.to_csv(pathname, index=False)
    return pathname


if __name__ == '__main__':
    # Example
    bold_data_01, bold_data_df_01 = load_data(
        "/Users/zhenrujia/Downloads/Sub_01_LSD_01.nii")
    bold_data_03, bold_data_df_03 = load_data(
        "/Users/zhenrujia/Downloads/Sub_01_LSD_03.nii")

    time_series_01, labels = parcellation("harvox", bold_data_01)
    time_series_03, labels = parcellation("harvox", bold_data_03)
    combined_series = combine_series(time_series_01, time_series_03)

    pathname = save_data(combined_series, "01", 1)
