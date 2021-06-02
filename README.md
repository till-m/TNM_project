# TNM_project

## Building the models:
1. Install the necessary python packages by running

        pip install -r requirements.txt


2. (Optional) Install the matlab engine by following the steps outlined [here](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html) **OR** execute the matlab script manually when prompted (not recommended and very tedious).

3. Download the datasets [ds003059](https://openneuro.org/datasets/ds003059/versions/1.0.0) and [ds00030](https://openneuro.org/datasets/ds000030/versions/1.0.0) to the `data/` directory.

4. (Optional) change the parcellation scheme by modifying the corresponding entry in `config.yml`. Options include:
    - `yeo` for [Yeo-2011](doi.org/10.1152/jn.00338.2011)
    - `harvox` for [Harvard-Oxford](https://neurovault.org/collections/262/)
    - `aal` for [Automated Anatomical Labeling](https://doi.org/10.1006%2Fnimg.2001.0978)
    -  `schaefer` for [Schafer-2018](doi.org/10.1093/cercor/bhx179)

5. Execute `main.py`.

## Reproducing the data analysis:
1. Make certain the necessary files are present: Either build the models yourself (follow the steps above) or download our pre-built models by running

        curl https://polybox.ethz.ch/index.php/s/jioOccqPFjU6cpk/download -o output_DCM.7z
    and unpacking to a folder called `output_DCM/`.

2. Run `analysis.m` using MATLAB.