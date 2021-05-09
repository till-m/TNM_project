% fix the random number generator
rng(2406,'twister') % is this still necessary?


% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));

% UNUSED UNUSED UNUSED
% specify the options for the rDCM analysis
options.SNR             = 3;
options.y_dt            = 0.5;
options.p0_all          = 0.15;  % single p0 value (for computational efficiency)
options.iter            = 100;
options.filter_str      = 5;
options.restrictInputs  = 1;
% UNUSED UNUSED UNUSED

fprintf('Load data\n')
DCM = load('.temp/in.mat');

disp(DCM);
% get time
currentTimer = tic;
DCM = tapas_rdcm_model_specification(DCM.Y, DCM.U, []);

disp(DCM)

type = 'r';
methods = 2;
[output, options] = tapas_rdcm_estimate(DCM, type, methods=methods);

% output elapsed time
toc(currentTimer)



%% visualize the results

% regions for which to plot traces
plot_regions = [1 12];

% visualize the results
tapas_rdcm_visualize(output, DCM, options, plot_regions, 1)



fprintf('Save to .temp/')
save('.temp/out.mat', 'output')

