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

%% Export results.
rDCM_output.meta = meta;

folder_path = "output_DCM/" + meta.scheme;
mkdir(folder_path);
fprintf("Saving to " + folder_path);
save(folder_path + "/" + meta.name, 'rDCM_output')
