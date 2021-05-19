% fix the random number generator
rng(2406,'twister')

% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));

fprintf('Load data\n')
in_mat = load('.temp/in.mat');

meta=in_mat.meta;
Y=in_mat.Y;
meta.regions = Y.name;

% get time
currentTimer = tic;
DCM = tapas_rdcm_model_specification(Y, [], []);

%disp(DCM)

type = 'r';
methods = 1;
[rDCM_output, options] = tapas_rdcm_estimate(DCM, type, [], methods);

% output elapsed time
toc(currentTimer)



%% visualize the results

% regions for which to plot traces
%plot_regions = [1 12];

% visualize the results
%tapas_rdcm_visualize(rDCM_output, DCM, options, plot_regions, 1)


%% Export results.
%A = rDCM_output.Ep.A;
A = rDCM_output;

folder_path = "output_DCM/" + meta.scheme;
mkdir(folder_path);
fprintf("Saving to " + folder_path);
%save(folder_path + "/" + meta.name, 'meta', 'A')
save(folder_path + "/" + meta.name, 'rDCM_output')

