% fix the random number generator
rng(2406,'twister') % is this still necessary?


% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));


fprintf('Load data\n')
Y = load('.temp/in.mat');

% get time
currentTimer = tic;
DCM = tapas_rdcm_model_specification(Y, [], []);

disp(DCM)

type = 'r';
methods = 1;
[output, options] = tapas_rdcm_estimate(DCM, type, methods=methods);

% output elapsed time
toc(currentTimer)



%% visualize the results

% regions for which to plot traces
%plot_regions = [1 12];

% visualize the results
%tapas_rdcm_visualize(output, DCM, options, plot_regions, 1)


%% Export results.
fprintf('Save to .temp/')
save('.temp/out.mat', 'output')

