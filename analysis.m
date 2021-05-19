clear all;
% fix the random number generator
rng(2406,'twister') % is this still necessary?

% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));
fprintf('Load data\n')

%% creating structs with all the rDCMoutputs from all the subjects: LSD_subjects_all and LSD_subjects_all
directory = "output_DCM/harvox/";

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
    
    LSD_subjects_all(i).name = allFileNames_PLCB{i};
    PLCB_subjects_all(i).rDCM_output = load(directory + allFileNames_PLCB{i}).rDCM_output;  
end

%% Analysis for parcellation
parcelations_matrix_size = size(PLCB_subjects_all(1).rDCM_output.Ep.A);
parcelation_n = parcelations_matrix_size(1);
meanMatrix = zeros(parcelations_matrix_size);
for i = 1 : n_subjects
    meanMatrix = meanMatrix + (LSD_subjects_all(i).rDCM_output.Ep.A - PLCB_subjects_all(i).rDCM_output.Ep.A);
end

meanMatrix = meanMatrix / n_subjects;

%from Valeries extra analysis
colormap('parula')
imagesc(meanMatrix)
colorbar
title('Average Difference from Yeo','FontSize',14)
axis square
caxis_range = [min(min(meanMatrix)) max(max(meanMatrix))];
caxis(caxis_range)
%caxis([-0.04 0.04]) %yeo
%caxis([-0.01 0.01]) %harvox
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

% flag the high average
significant = meanMatrix;
%significant(abs(significant)<=0.02) = 1e-9; %yeo threshold
significant(abs(significant)<=0.01) = 1e-9; %harvox threshold
colormap('parula')
imagesc(significant)
colorbar
title('Average Difference','FontSize',14)
axis square
caxis(caxis_range)
%caxis([-0.04 0.04]) %yeo threshold
%caxis([-0.01 0.01]) %harvox threshold
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
set(gca,'xtick',[1:parcelation_n])
set(gca,'ytick',[1:parcelation_n])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

%% std deviation stuff
% calculate the std
temp = 0;
for i = 1 : n_subjects
    A_difference = (LSD_subjects_all(i).rDCM_output.Ep.A - PLCB_subjects_all(i).rDCM_output.Ep.A);
    temp = temp + (A_difference - meanMatrix).^2;
end

std = sqrt(temp./n_subjects);
std(abs(std)<=0.0126) = 1e-9;
colormap('parula')
imagesc(std)
colorbar

% flag the high average
significant = meanMatrix;
significant(abs(significant)<=0.007) = 1e-9;
colormap('parula')
imagesc(significant)
colorbar
title('Average Difference','FontSize',14)
axis square
%caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
caxis([-0.01 0.01]) 
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
%set(gca,'xtick',[1:48])
%set(gca,'ytick',[1:48])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

%% Significance test
temp = zeros(1,n_subjects);
lst = zeros(1, n_subjects);
for i = 1:parcelation_n
    for j = 1:parcelation_n

        %lst = [Sub01_diff_harvox(i,j) Sub02_diff_harvox(i,j) Sub03_diff_harvox(i,j) Sub04_diff_harvox(i,j) Sub06_diff_harvox(i,j) Sub09_diff_harvox(i,j) Sub10_diff_harvox(i,j) Sub11_diff_harvox(i,j) Sub12_diff_harvox(i,j) Sub13_diff_harvox(i,j) Sub15_diff_harvox(i,j) Sub17_diff_harvox(i,j) Sub18_diff_harvox(i,j) Sub19_diff_harvox(i,j) Sub20_diff_harvox(i,j)];
        for k = 1:n_subjects
            A_difference = (LSD_subjects_all(k).rDCM_output.Ep.A - PLCB_subjects_all(k).rDCM_output.Ep.A);
            lst(1,k) = A_difference(i,j);
        end
        %temp = cat(1,temp,lst);
        temp = cat(1,temp,lst);
    end
end
final = temp(2:end,:);

% find the indices of 20 flagged regions
linearIndexes = find(abs(significant) > 1e-9); 
[rows, columns] = ind2sub(size(significant), linearIndexes);


% 20 vs all
h_lst = zeros(parcelation_n);
for i = 1:20
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
    h_lst(r_ind,columns(i)) = mode(temp);
end


% all vs all
%h_lst = zeros(48);
%for i = 1:48*48
%    temp = zeros(1,48*48);
%    for j = 1:48*48
%        [h,p] = ttest2(final(i,:), final(j,:));
%        temp(j) = h;
        %p_lst(i,j+1) = p;
%    end
%    r_ind = ceil(i/48);
%    if rem(i,48) == 0
%        c_ind = 48;
%    else
%        c_ind = rem(i,48);
%    end
%    h_lst(r_ind,c_ind) = mode(temp);
%end


imagesc(h_lst) 
map = [0.2 0.1 0.5
   0 1 1];
colormap(map)
%set(gca,'xtick',[1:48])
%set(gca,'ytick',[1:48])
title('Significance test for the mean','FontSize',14)
axis image; hold on
for i=1:20
    plot(columns(i),rows(i),'marker','x','LineWidth',3,'Color','red')
    hold on
end
L = line(ones(2), ones(2));
set(L, {'Color'}, num2cell(map, 2))
legend(L, {'0','1'})
hold off


