clear all;
% fix the random number generator
rng(2406,'twister')

% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));
fprintf('Load data\n')

% Basic settings necessary to set before start
mode = 0; %1: SCZ, 0:LSD
%directory = "output_DCM/SCZ/yeo/";
directory = "yeo/";
signficance_threshold = 1e-3; %yeo threshold 0.02, harvox threshold 0.01
% End Basic settings necessary to set before start

if mode == 0
    
    %% creating structs with all the rDCMoutputs from all the subjects: LSD_subjects_all and LSD_subjects_all
    listing_LSD = dir(directory+'*LSD.mat');
    listing_PLCB = dir(directory+'*PLCB.mat');
    allFileNames_LSD = {listing_LSD(:).name};
    allFileNames_PLCB = {listing_PLCB(:).name};
    n_subjects = length(allFileNames_LSD);

    for k = 1 : n_subjects
      fprintf('allFileNames_LSD{%d} = %s\n', k, allFileNames_LSD{k});
      fprintf('allFileNames_PLCB{%d} = %s\n', k, allFileNames_PLCB{k});
    end

    for i = 1 : n_subjects
        LSD_subjects_all(i).name = allFileNames_LSD{i};
        LSD_subjects_all(i).rDCM_output = load(directory + allFileNames_LSD{i}).rDCM_output;

        PLCB_subjects_all(i).name = allFileNames_PLCB{i};
        PLCB_subjects_all(i).rDCM_output = load(directory + allFileNames_PLCB{i}).rDCM_output;  
    end

    %% Analysis for parcellation
    all_regions = cellstr(PLCB_subjects_all(1).rDCM_output.meta.regions);
    
    parcelations_matrix_size = size(PLCB_subjects_all(1).rDCM_output.Ep.A);
    parcelation_n = parcelations_matrix_size(1);
    meanMatrix = zeros(parcelations_matrix_size);

    for i = 1 : n_subjects
        meanMatrix = meanMatrix + (LSD_subjects_all(i).rDCM_output.Ep.A - PLCB_subjects_all(i).rDCM_output.Ep.A);
        diff = LSD_subjects_all(i).rDCM_output.Ep.A - PLCB_subjects_all(i).rDCM_output.Ep.A;
        fileName = sprintf('diff_LSD_%d.mat', i);
        save(fileName,'diff')
    end

    meanMatrix = meanMatrix ./ n_subjects;
    
    std = 0;
    for i = 1 : n_subjects
        A_difference = (LSD_subjects_all(i).rDCM_output.Ep.A - PLCB_subjects_all(i).rDCM_output.Ep.A);
        std = std + (A_difference - meanMatrix).^2;
    end
    std = sqrt(std./n_subjects);
    
    % create a matrix of all_entries_of_matrix_A * n_subjects for the
    % significance test
    temp = zeros(1,n_subjects);
    for i = 1:parcelation_n
        for j = 1:parcelation_n
            lst = zeros(1,n_subjects);
            for k = 1:n_subjects
                fileName = sprintf('diff_LSD_%d.mat', k);
                diffMatrix = load(fileName).diff;
                lst(k) = diffMatrix(i,j);
            end
            temp = cat(1,temp,lst);
        end
    end
    final = temp(2:end,:);


else
    listing_SCZ = dir(directory+'*SCZ.mat');
    listing_CTRL = dir(directory+'*CTRL.mat');
    allFileNames_SCZ = {listing_SCZ(:).name};
    allFileNames_CTRL = {listing_CTRL(:).name};
    n_subjects_SCZ = length(allFileNames_SCZ);    
    n_subjects_CTRL = length(allFileNames_CTRL); 


    parcelations_matrix_size = size(load(directory + allFileNames_SCZ{1}).rDCM_output.Ep.A);
    parcelation_n = parcelations_matrix_size(1);
    SCZ_meanMatrix = zeros(parcelations_matrix_size);
    CTRL_meanMatrix = zeros(parcelations_matrix_size); 

    for i = 1 : n_subjects_SCZ
        SCZ_subjects_all(i).name = allFileNames_SCZ{i};
        SCZ_subjects_all(i).rDCM_output = load(directory + allFileNames_SCZ{i}).rDCM_output;
        SCZ_meanMatrix = SCZ_meanMatrix + SCZ_subjects_all(i).rDCM_output.Ep.A;
    end
    for i = 1 : n_subjects_CTRL    
        CTRL_subjects_all(i).name = allFileNames_CTRL{i};
        CTRL_subjects_all(i).rDCM_output = load(directory + allFileNames_CTRL{i}).rDCM_output;
        CTRL_meanMatrix = CTRL_meanMatrix + CTRL_subjects_all(i).rDCM_output.Ep.A;
    end

    all_regions = cellstr(CTRL_subjects_all(1).rDCM_output.meta.regions);

    SCZ_meanMatrix = SCZ_meanMatrix ./ n_subjects_SCZ;
    CTRL_meanMatrix = CTRL_meanMatrix ./ n_subjects_CTRL;
    %meanMatrix = SCZ_meanMatrix - CTRL_meanMatrix;

    meanMatrix = zeros(parcelations_matrix_size);
    for i = 1 : n_subjects_SCZ
        SCZ_subjects_all(i).name = allFileNames_SCZ{i};
        SCZ_subjects_all(i).rDCM_output = load(directory + allFileNames_SCZ{i}).rDCM_output;
        meanMatrix = meanMatrix + (SCZ_subjects_all(i).rDCM_output.Ep.A - CTRL_meanMatrix);
        %diff = SCZ_subjects_all(i).rDCM_output.Ep.A - CTRL_meanMatrix;
        %fileName = sprintf('diff_SCZ_%d.mat', i);
        %save(fileName,'diff')
    end
    meanMatrix = meanMatrix ./ n_subjects_SCZ;
    
    
    % create a matrix of all_entries_of_matrix_A * n_subjects for the
    % significance test
    temp = zeros(1,n_subjects_SCZ);
    for i = 1:parcelation_n
        for j = 1:parcelation_n
            lst = zeros(1,n_subjects_SCZ);
            for k = 1:n_subjects_SCZ
                fileName = sprintf('diff_SCZ_%d.mat', k);
                diffMatrix = load(fileName).diff;
                lst(k) = diffMatrix(i,j);
            end
            temp = cat(1,temp,lst);
        end
    end
    final = temp(2:end,:);

end


regions_to_analyse = all_regions; %{'Planum Temporale'; 'Supracalcarine Cortex'; 'Occipital Pole'}
meanMatrix = select_regions_by_name(regions_to_analyse, all_regions, meanMatrix);
%from Valeries extra analysis
colormap('parula')
imagesc(meanMatrix)
colorbar
title('Average Difference','FontSize',14)
axis square
caxis_range = [min(min(meanMatrix)) max(max(meanMatrix))];
caxis(caxis_range)
%caxis([-0.04 0.04]) %yeo LSD
%caxis([-0.01 0.01]) %harvox LSD
%caxis([-0.025 0.025]) %yeo SCZ
caxis([-0.008 0.008]) %harvox SCZ
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

% flag the high average
significant = meanMatrix;
%significant(abs(significant)<=signficance_threshold) = 1e-9;
%significant(abs(significant)<=0.02) = 1e-9; %yeo threshold LSD
%significant(abs(significant)<=0.007) = 1e-9; %harvox threshold LSD
%significant(abs(significant)<=0.015) = 1e-9; %yeo threshold SCZ
significant(abs(significant)<=0.0055) = 1e-9; %harvox threshold LSD
colormap('parula')
imagesc(significant)
colorbar
title('Average Difference','FontSize',14)
axis square
%caxis(caxis_range)
%caxis([-0.04 0.04]) %yeo threshold LSD
%caxis([-0.01 0.01]) %harvox threshold LSD
caxis([-0.025 0.025]) %yeo SCZ
%caxis([-0.008 0.008]) %harvox SCZ
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
set(gca,'xtick',[1:parcelation_n])
set(gca,'ytick',[1:parcelation_n])
%set(gca,'xticklabels', regions_to_analyse)
%set(gca,'yticklabels', regions_to_analyse)
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)
 

% find the indices of flagged regions
linearIndexes = find(abs(significant) > 1e-9); 
[rows, columns] = ind2sub(size(significant), linearIndexes);

n_significant = size(rows);
n_significant = n_significant(1);

% significant vs all
h_lst = zeros(parcelation_n);
for i = 1:n_significant
    r_ind = rows(i);
    temp = zeros(1,parcelation_n*parcelation_n);
    for j = 1:parcelation_n*parcelation_n
        [h,p] = ttest2(final((r_ind-1)*parcelation_n + columns(i),:), final(j,:));
        if (r_ind-1)*parcelation_n + columns(i) == j
            temp(j) = 100;
        else
            temp(j) = h;
        end
        %p_lst(i,j+1) = p;
    end
    %h_lst(r_ind,columns(i)) = mode(temp);
    %disp((sum(temp)-100)/(parcelation_n*parcelation_n-1))
    %if (sum(temp)-100)/(parcelation_n*parcelation_n-1) > 0.6 % for LSD
    if (sum(temp)-100)/(parcelation_n*parcelation_n-1) > 0.75 % for SCZ
        h_lst(r_ind,columns(i)) = 1;
    else 
        h_lst(r_ind,columns(i)) = 0;
    end
end


imagesc(h_lst) 
map = [0.2 0.1 0.5
   0 1 1];
colormap(map)
%set(gca,'xtick',[1:48])
%set(gca,'ytick',[1:48])
title('Significance test for the mean','FontSize',14)
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)
axis image; hold on
for i=1:n_significant
    plot(columns(i),rows(i),'marker','x','LineWidth',3,'Color','red')
    hold on
end
L = line(ones(2), ones(2));
set(L, {'Color'}, num2cell(map, 2))
legend(L, {'0','1'})
hold off




% 
% %% std deviation stuff
% % calculate the std
% temp = 0;
% for i = 1 : n_subjects
%     A_difference = (LSD_subjects_all(i).rDCM_output.Ep.A - PLCB_subjects_all(i).rDCM_output.Ep.A);
%     temp = temp + (A_difference - meanMatrix).^2;
% end
% 
% std = sqrt(temp./n_subjects);
% std(abs(std)<=0.0126) = 1e-9;
% colormap('parula')
% imagesc(std)
% colorbar
% 
% % flag the high average
% significant = meanMatrix;
% significant(abs(significant)<=0.007) = 1e-9;
% colormap('parula')
% imagesc(significant)
% colorbar
% title('Average Difference','FontSize',14)
% axis square
% %caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
% caxis([-0.01 0.01]) 
% %set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
% %set(gca,'xtick',[1:48])
% %set(gca,'ytick',[1:48])
% %set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
% xlabel('region (from)','FontSize',12)
% ylabel('region (to)','FontSize',12)
% 
% %% Significance test
% temp = zeros(1,n_subjects);
% lst = zeros(1, n_subjects);
% for i = 1:parcelation_n
%     for j = 1:parcelation_n
% 
%         %lst = [Sub01_diff_harvox(i,j) Sub02_diff_harvox(i,j) Sub03_diff_harvox(i,j) Sub04_diff_harvox(i,j) Sub06_diff_harvox(i,j) Sub09_diff_harvox(i,j) Sub10_diff_harvox(i,j) Sub11_diff_harvox(i,j) Sub12_diff_harvox(i,j) Sub13_diff_harvox(i,j) Sub15_diff_harvox(i,j) Sub17_diff_harvox(i,j) Sub18_diff_harvox(i,j) Sub19_diff_harvox(i,j) Sub20_diff_harvox(i,j)];
%         for k = 1:n_subjects
%             A_difference = (LSD_subjects_all(k).rDCM_output.Ep.A - PLCB_subjects_all(k).rDCM_output.Ep.A);
%             lst(1,k) = A_difference(i,j);
%         end
%         %temp = cat(1,temp,lst);
%         temp = cat(1,temp,lst);
%     end
% end
% final = temp(2:end,:);
% 
% % find the indices of 20 flagged regions
% linearIndexes = find(abs(significant) > 1e-9); 
% [rows, columns] = ind2sub(size(significant), linearIndexes);
% 
% 
% % 20 vs all
% h_lst = zeros(parcelation_n);
% for i = 1:20
%     r_ind = rows(i);
%     temp = zeros(1,parcelation_n*parcelation_n);
%     for j = 1:parcelation_n*parcelation_n
%         [h,p] = ttest2(final((r_ind-1)*parcelation_n + columns(i),:), final(j,:));
%         if (r_ind-1)*parcelation_n + columns(i) == j
%             temp(j) = 100;
%         else
%             temp(j) = h;
%         end
%         %p_lst(i,j+1) = p;
%     end
%     h_lst(r_ind,columns(i)) = mode(temp);
% end
% 
% 
% % all vs all
% %h_lst = zeros(48);
% %for i = 1:48*48
% %    temp = zeros(1,48*48);
% %    for j = 1:48*48
% %        [h,p] = ttest2(final(i,:), final(j,:));
% %        temp(j) = h;
%         %p_lst(i,j+1) = p;
% %    end
% %    r_ind = ceil(i/48);
% %    if rem(i,48) == 0
% %        c_ind = 48;
% %    else
% %        c_ind = rem(i,48);
% %    end
% %    h_lst(r_ind,c_ind) = mode(temp);
% %end
% 
% 
% imagesc(h_lst) 
% map = [0.2 0.1 0.5
%    0 1 1];
% colormap(map)
% %set(gca,'xtick',[1:48])
% %set(gca,'ytick',[1:48])
% title('Significance test for the mean','FontSize',14)
% axis image; hold on
% for i=1:20
%     plot(columns(i),rows(i),'marker','x','LineWidth',3,'Color','red')
%     hold on
% end
% L = line(ones(2), ones(2));
% set(L, {'Color'}, num2cell(map, 2))
% legend(L, {'0','1'})
% hold off
% 
% 

%% clustering of SCZ data
out = cell(1,n_subjects_SCZ);
visual = cell(1,n_subjects_SCZ);
num = zeros(1,n_subjects_SCZ);
for k = 1:n_subjects_SCZ
    fileName = sprintf('yeo/diff_SCZ_%d.mat', k);
    out{k} = load(fileName).diff;
    visual{k} = out{k}(:,1:2);
    visual{k}(abs(visual{k})<=0.0162) = 1e-9; 
    linearIndexes = find(abs(visual{k}) > 1e-9); 
    num(k) = size(linearIndexes,1);
end

plot(num,'-o')
yline(21,'r--', 'LineWidth', 3)
xlabel("Index of Schizophrenia Patients")
ylabel("Number of entries with values higher than the threshold")



visual_lst = find(num > 20);
absMeanMatrix = zeros(parcelation_n);
for i= 1:size(visual_lst,2)
    absMeanMatrix = absMeanMatrix + abs(out{visual_lst(i)});
end
absMeanMatrix = absMeanMatrix ./ size(visual_lst,2);
            

temp = zeros(1,size(visual_lst,2));
for i = 1:parcelation_n
    for j = 1:parcelation_n
        lst = zeros(1,size(visual_lst,2));
        for k = 1:size(visual_lst,2)
            ind = visual_lst(k);
            lst(k) = abs(out{ind}(i,j));
        end
        temp = cat(1,temp,lst);
    end
end
final_visual = temp(2:end,:);


colormap('parula')
imagesc(absMeanMatrix)
colorbar
caxis([0 0.05]) 
title('Average Absolute Difference for Subjects with High Visual Cortex Impacts','FontSize',12)
axis square
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)


% flag the high average
significant = absMeanMatrix;
significant(significant<=0.036) = 1e-9; 
colormap('parula')
imagesc(significant)
colorbar
title('Average Absolute Difference for Subjects with High Visual Cortex Impacts','FontSize',12)
axis square
caxis([0 0.05]) 
set(gca,'xtick',[1:parcelation_n])
set(gca,'ytick',[1:parcelation_n])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)
 

full_lst = [1:n_subjects_SCZ];
notVisual_lst = setdiff(full_lst, visual_lst); 

absMeanMatrix_notVisual = zeros(parcelation_n);
for i= 1:size(notVisual_lst,2)
    absMeanMatrix_notVisual = absMeanMatrix_notVisual + abs(out{notVisual_lst(i)});
end
absMeanMatrix_notVisual = absMeanMatrix_notVisual ./ size(notVisual_lst,2);

temp = zeros(1,size(notVisual_lst,2));
for i = 1:parcelation_n
    for j = 1:parcelation_n
        lst = zeros(1,size(notVisual_lst,2));
        for k = 1:size(notVisual_lst,2)
            ind = notVisual_lst(k);
            lst(k) = abs(out{ind}(i,j));
        end
        temp = cat(1,temp,lst);
    end
end
final_notVisual = temp(2:end,:);


colormap('parula')
imagesc(absMeanMatrix_notVisual)
colorbar
caxis([0 0.05]) 
title('Average Absolute Difference for Subjects with Low Visual Cortex Impacts','FontSize',12)
axis square
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)



% flag the high average
significant_notVisual = absMeanMatrix_notVisual;
significant_notVisual(significant_notVisual<=0.036) = 1e-9; 
colormap('parula')
imagesc(significant_notVisual)
colorbar
title('Average Absolute Difference for Subjects with Low Visual Cortex Impacts','FontSize',12)
axis square
caxis([0 0.05]) 
set(gca,'xtick',[1:parcelation_n])
set(gca,'ytick',[1:parcelation_n])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)



%% functions to select regions
function res = select_regions_by_number(regions_idx, matrix)
    res = matrix(regions_idx,regions_idx);
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


   