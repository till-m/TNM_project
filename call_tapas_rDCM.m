% fix the random number generator
rng(2406,'twister') % is this still necessary?


% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));

% specify the options for the rDCM analysis
options.SNR             = 3;
options.y_dt            = 0.5;
options.p0_all          = 0.15;  % single p0 value (for computational efficiency)
options.iter            = 100;
options.filter_str      = 5;
options.restrictInputs  = 1;


fprintf('Load data\n')
DCM = '.temp/in.mat'; %load('.temp/in.mat').DCM;

% get time
currentTimer = tic;

type = 'r';
methods = 1;
out = tapas_rdcm_estimate(DCM, type, options, methods);
fprintf('Save to .temp')
save('.temp/rDCM.mat', 'out')
