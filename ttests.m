% does analysis

% fix the random number generator
rng(2406,'twister')


%% MAIN

analysis('yeo', [-0.025 0.025], 1)
analysis('schaefer', [-0.002 0.002], 1)

function analysis(name, caxis_range, FDR_correction)
    LSD_subjects = load_data("output_DCM/" +name +"/", "LSD");
    PLCB_subjects = load_data("output_DCM/" +name +"/", "PLCB");
    SCZ_subjects = load_data("output_DCM/" +name +"/", "SCZ");
    CTRL_subjects = load_data("output_DCM/" +name +"/", "CTRL");
    ticklabels = cellstr(LSD_subjects(1).rDCM_output.meta.regions);
    ttest_wrapper(LSD_subjects, PLCB_subjects, SCZ_subjects, CTRL_subjects, FDR_correction, caxis_range);%, ticklabels)
end


%% auxiliary function definitions
function ttest_wrapper(LSD_subjects, PLCB_subjects, SCZ_subjects, CTRL_subjects, FDR_correction, caxis_range, ticklabels)
    % bad way of making things optional
    if ~(exist('ticklabels', 'var'))
        ticklabels = [];
    end

    % t-test LSD vs. PLCB
    tt = t_test(LSD_subjects, PLCB_subjects, FDR_correction);
    plot_significance(tt,'Significance LSD/PLCB', ticklabels);
    
    % t-test SCZ vs. CTRL
    tt = t_test(SCZ_subjects, CTRL_subjects, FDR_correction);
    plot_significance(tt,'Significance SCZ/CTRL', ticklabels);
    
    % t-test SCZ vs. CTRL
    diff1 = unpaired_diff(LSD_subjects, PLCB_subjects);
    diff2 = unpaired_diff(SCZ_subjects, CTRL_subjects);
    
    tt = t_test(diff1, diff2, FDR_correction);
    plot_significance(tt,'Significance LSD-PLCB_{avg}/SCZ-CTRL_{avg}', ticklabels);
    
    [cor,Pval] = correlation(diff1,diff2);
    disp(cor)
    disp(Pval)
    
    
    %LSD_SCZ = mean(concat_subjects(diff1).') - mean(concat_subjects(diff2).');
    %shape = size(LSD_subjects(1).rDCM_output.Ep.A);
    %LSD_SCZ = reshape(LSD_SCZ,shape);
    %LSD_SCZ(tt==0) = 0; 
    %plot_matrix_diff(LSD_SCZ,'Average Differences [LSD-PLCB_{avg}] - [SCZ-CTRL_{avg}]', [], ticklabels);


    avg_and_plot_matrix(diff1, tt, 'Average LSD-PLCB_{avg}', caxis_range, ticklabels);
    avg_and_plot_matrix(diff2, tt, 'Average SCZ-CTRL_{avg}', caxis_range, ticklabels);
    
end

function res = t_test(subjects1, subjects2, FDR_correction)
    shape = size(subjects1(1).rDCM_output.Ep.A);
    [res,p] = ttest2(concat_subjects(subjects1).', concat_subjects(subjects2).'); 
    % FDR correction
    if FDR_correction == 1
        [~,q] = mafdr(p);
        res = q <= 0.05;
    end
    res = reshape(res, shape);
end

function res = concat_subjects(subjects)
    n_subjects = size(subjects, 2);
    res = [];
    for i = 1:n_subjects
        col = subjects(i).rDCM_output.Ep.A(:);
        res = [res, col];
    end
end

function diff = paired_diff(subjects1, subjects2)
    n_subjects = size(subjects1, 2);
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
    for i = 2:n_subjects
        res = res + subjects(i).rDCM_output.Ep.A;
    end
    res = res ./ n_subjects;
end

function [res,P] = correlation(subjects1,subjects2)
    size(mean(concat_subjects(subjects1).'))
    [res,P] = corrcoef(mean(concat_subjects(subjects1).'), mean(concat_subjects(subjects2).'));
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

function plot_significance(matrix, plot_title, ticklabels)
    figure()

    map = [0.2 0.1 0.5
           0 1 1];
    colormap(map)
    imagesc(matrix)
    %colorbar
    title(plot_title, 'FontSize', 14)
    axis square

    xlabel('region (from)','FontSize',12)
    ylabel('region (to)','FontSize',12)
    %set(gca,'xtick',[1:size(matrix,1)])
    %set(gca,'ytick',[1:size(matrix,1)])
    L = line(ones(2), ones(2));
    set(L, {'Color'}, num2cell(map, 2))
    legend(L, {'Not significant','Significant'},'Location','northeastoutside')
    if ~(size(ticklabels,1)==0)
        set(gca,'xtick',1:size(matrix,1))
        set(gca,'ytick',1:size(matrix,1))
        set(gca,'xticklabels', ticklabels)
        set(gca,'yticklabels', ticklabels)
    end
    shg
end

function avg_and_plot_matrix(diff, tt_result, plot_title, caxis_range, ticklabels)
    avg = mean(concat_subjects(diff).');
    shape = size(tt_result);
    avg = reshape(avg,shape);
    avg(tt_result==0) = 0;

    figure()

    colormap('parula')
    imagesc(avg)
    colorbar

    title(plot_title, 'FontSize', 14)
    axis square
    if ~(size(caxis_range,1)==0)
        caxis(caxis_range)
    end
    xlabel('region (from)','FontSize',12)
    ylabel('region (to)','FontSize',12)

    if ~(size(ticklabels,1)==0)
        set(gca,'xtick',1:size(avg,1))
        set(gca,'ytick',1:size(avg,1))
        set(gca,'xticklabels', ticklabels)
        set(gca,'yticklabels', ticklabels)
    end
    shg
end