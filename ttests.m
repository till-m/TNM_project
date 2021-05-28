% does analysis

% fix the random number generator
rng(2406,'twister')

%% MAIN

yeo()
schaefer()



%% yeo analysis
function yeo()
    LSD_subjects = load_data("output_DCM/yeo/", "LSD");
    PLCB_subjects = load_data("output_DCM/yeo/", "PLCB");
    SCZ_subjects = load_data("output_DCM/yeo/", "SCZ");
    CTRL_subjects = load_data("output_DCM/yeo/", "CTRL");
    ttest_wrapper(LSD_subjects, PLCB_subjects, SCZ_subjects, CTRL_subjects, 0)
    
end


%% schaefer analysis
function schaefer()
    LSD_subjects = load_data("output_DCM/schaefer/", "LSD");
    PLCB_subjects = load_data("output_DCM/schaefer/", "PLCB");
    SCZ_subjects = load_data("output_DCM/schaefer/", "SCZ");
    CTRL_subjects = load_data("output_DCM/schaefer/", "CTRL");
    ttest_wrapper(LSD_subjects, PLCB_subjects, SCZ_subjects, CTRL_subjects, 0)
end

%% function definitions
function ttest_wrapper(LSD_subjects, PLCB_subjects, SCZ_subjects, CTRL_subjects, ticklabels)
    % t-test LSD vs. PLCB
    tt = t_test(LSD_subjects, PLCB_subjects);
    plot_matrix(tt,'Significance LSD/PLCB', 0, ticklabels);
    
    % t-test SCZ vs. CTRL
    tt = t_test(SCZ_subjects, CTRL_subjects);
    plot_matrix(tt,'Significance LSD/PLCB', 0, ticklabels);
    
    % t-test SCZ vs. CTRL
    diff1 = unpaired_diff(LSD_subjects, PLCB_subjects);
    diff2 = unpaired_diff(SCZ_subjects, CTRL_subjects);
    
    tt = t_test(diff1, diff2);
    plot_matrix(tt,'Significance LSD-PLCB_{avg}/SCZ-CTRL_{avg}', 0, ticklabels);
end

function res = t_test(subjects1, subjects2)
    shape = size(subjects1(1).rDCM_output.Ep.A);
    res = ttest2(concat_subjects(subjects1).', concat_subjects(subjects2).'); 
    res = reshape(res, shape);
end

function res = concat_subjects(subjects)
    n_subjects = size(subjects, 2);
    res = [];
    for i = 1:n_subjects
        disp(subjects(i).name)
        col = subjects(i).rDCM_output.Ep.A(:);
        res = [res, col];
    end
end

function diff = paired_diff(subjects1, subjects2)
    n_subjects = size(subjects, 2);
    for i = 1:n_subjects
        diff(i).name = subjects1(i).name;
        diff(i).rDCM_output.Ep.A = subjects1(i).rDCM_output.Ep.A - subjects2(i).rDCM_output.Ep.A;
    end
end

function diff = unpaired_diff(subjects1, subjects2)
    n_subjects1 = size(subjects1, 2);
    subjects2_avg = average_over_subjects(subjects2);
    for i = 1:n_subjects1
        diff(i).name = subjects1(i).name;
        diff(i).rDCM_output.Ep.A = subjects1(i).rDCM_output.Ep.A - subjects2_avg;
    end
end

function res = average_over_subjects(subjects)
    n_subjects = size(subjects, 2);
    res = subjects(1).rDCM_output.Ep.A;
    for i = 1:n_subjects
        res = res + subjects(i).rDCM_output.Ep.A;
    end
    res = res ./ n_subjects;
end

function all_subjects = load_data(directory, type)
    files = dir(directory + "*" + type + ".mat");
    all_file_names = {files(:).name};
    n_subjects = length(all_file_names);
    all_subjects = [];
    for i = 1 : n_subjects
        all_subjects(i).name = all_file_names{i};
        all_subjects(i).rDCM_output = load(directory + all_file_names{i}).rDCM_output;
    end
end

function subjects = select_regions_by_name_multiple(regions, subjects)
    all_regions = cellstr(subjects(1).rDCM_output.meta.regions);
    n_subjects = size(subjects, 2);
    for i = 1:n_subjects
        subjects(i).rDCM_output.Ep.A = select_regions_by_name(regions, all_regions, subjects(i).rDCM_output.Ep.A);
    end
end

function res = select_regions_by_name(regions, all_regions, matrix)
    n_regions = size(regions, 1);
    idx = [];
    for i = 1:n_regions
        region = regions(i);
        X = contains(all_regions, region);
        idx(end+1) = find(X);
    end
    res = matrix(idx,idx);
end

function plot_matrix(matrix, plot_title, caxis_range, ticklabels)
    figure()
    % no keyword arguments, what a joke of a language.
    % set caxis_range and ticklabels to 0 if NA
    colormap('parula')
    imagesc(matrix)
    colorbar
    title(plot_title, 'FontSize', 14)
    axis square
    if ~(caxis_range==0)
        caxis(caxis_range)
    end
    set(gca,'xtick',1:size(matrix,1))
    set(gca,'ytick',1:size(matrix,1))
    xlabel('region (from)','FontSize',12)
    ylabel('region (to)','FontSize',12)
    if ~(ticklabels==0)
        set(gca,'xticklabels', ticklabels)
        set(gca,'yticklabels', ticklabels)
    end
    shg
end