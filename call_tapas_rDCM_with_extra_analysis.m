% fix the random number generator
rng(2406,'twister') % is this still necessary?

% get path of rDCM toolbox
P        = mfilename('fullpath');
rDCM_ind = strfind(P,fullfile('rDCM','code'));

fprintf('Load data\n')
%Y = load('.temp/in.mat');
Y = load('Y.mat');
disp(Y.y)

% get time
currentTimer = tic;
DCM = tapas_rdcm_model_specification(Y, [], []);
%disp(DCM)
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


colormap('parula')
imagesc(output.Ep.A)
colorbar
title('estimated','FontSize',14)
axis square
%caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
caxis([-0.1 0.1]) 
set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

%% Export results.
fprintf('Save to .temp/')
save('.temp/out.mat', 'output')


%% Obtain results for difference.
A = output.Ep.A;
save('A.mat','A')

PLCB_A = load('A.mat');
% difference = LSD effect - PLCB effect
diff = output.Ep.A - PLCB_A.A;
save('Sub20_diff_harvox.mat','diff')


colormap('parula')
imagesc(diff)
colorbar
title('Difference','FontSize',14)
axis square
%caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
caxis([-0.04 0.04]) 
set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)



%% Analysis for the Yeo parcellation
Sub01_diff = load('Sub01_diff.mat').diff;
Sub02_diff = load('Sub02_diff.mat').diff;
Sub03_diff = load('Sub03_diff.mat').diff;
Sub04_diff = load('Sub04_diff.mat').diff;
Sub06_diff = load('Sub06_diff.mat').diff;
Sub09_diff = load('Sub09_diff.mat').diff;
Sub10_diff = load('Sub10_diff.mat').diff;
Sub11_diff = load('Sub11_diff.mat').diff;
Sub12_diff = load('Sub12_diff.mat').diff;
Sub13_diff = load('Sub13_diff.mat').diff;
Sub15_diff = load('Sub15_diff.mat').diff;
Sub17_diff = load('Sub17_diff.mat').diff;
Sub18_diff = load('Sub18_diff.mat').diff;
Sub19_diff = load('Sub19_diff.mat').diff;
Sub20_diff = load('Sub20_diff.mat').diff;

meanMatrix = (Sub01_diff + Sub02_diff + Sub03_diff + Sub04_diff + Sub06_diff + Sub09_diff + Sub10_diff + Sub11_diff + Sub12_diff + Sub13_diff + Sub15_diff + Sub17_diff + Sub18_diff + Sub19_diff + Sub20_diff)/15;

colormap('parula')
imagesc(meanMatrix)
colorbar
title('Average Difference from Yeo','FontSize',14)
axis square
%caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
caxis([-0.04 0.04]) 
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

% flag the high average
significant = meanMatrix;
significant(abs(significant)<=0.02) = 1e-9;
colormap('parula')
imagesc(significant)
colorbar
title('Average Difference from Yeo','FontSize',14)
axis square
%caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
caxis([-0.04 0.04]) 
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
set(gca,'xtick',[1:17])
set(gca,'ytick',[1:17])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)


%% Analysis for the Harvox parcellation
Sub01_diff_harvox = load('Sub01_diff_harvox.mat').diff;
Sub02_diff_harvox = load('Sub02_diff_harvox.mat').diff;
Sub03_diff_harvox = load('Sub03_diff_harvox.mat').diff;
Sub04_diff_harvox = load('Sub04_diff_harvox.mat').diff;
Sub06_diff_harvox = load('Sub06_diff_harvox.mat').diff;
Sub09_diff_harvox = load('Sub09_diff_harvox.mat').diff;
Sub10_diff_harvox = load('Sub10_diff_harvox.mat').diff;
Sub11_diff_harvox = load('Sub11_diff_harvox.mat').diff;
Sub12_diff_harvox = load('Sub12_diff_harvox.mat').diff;
Sub13_diff_harvox = load('Sub13_diff_harvox.mat').diff;
Sub15_diff_harvox = load('Sub15_diff_harvox.mat').diff;
Sub17_diff_harvox = load('Sub17_diff_harvox.mat').diff;
Sub18_diff_harvox = load('Sub18_diff_harvox.mat').diff;
Sub19_diff_harvox = load('Sub19_diff_harvox.mat').diff;
Sub20_diff_harvox = load('Sub20_diff_harvox.mat').diff;

meanMatrix_harvox = (Sub01_diff_harvox + Sub02_diff_harvox + Sub03_diff_harvox + Sub04_diff_harvox + Sub06_diff_harvox + Sub09_diff_harvox + Sub10_diff_harvox + Sub11_diff_harvox + Sub12_diff_harvox + Sub13_diff_harvox + Sub15_diff_harvox + Sub17_diff_harvox + Sub18_diff_harvox + Sub19_diff_harvox + Sub20_diff_harvox)/15;
colormap('parula')
imagesc(meanMatrix_harvox)
colorbar
title('Average Difference from Harvard-Oxford','FontSize',14)
axis square
%caxis([-1*max(max(abs(output.Ep.A-diag(diag(output.Ep.A))))) max(max(abs(output.Ep.A-diag(diag(output.Ep.A)))))])
caxis([-0.01 0.01]) 
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)



% calculate the std
temp= (Sub01_diff_harvox - meanMatrix_harvox ) .^2 + (Sub02_diff_harvox - meanMatrix_harvox ) .^2 + (Sub03_diff_harvox - meanMatrix_harvox ) .^2 + (Sub04_diff_harvox - meanMatrix_harvox ) .^2;
temp = temp + (Sub06_diff_harvox - meanMatrix_harvox ) .^2 + (Sub09_diff_harvox- meanMatrix_harvox ) .^2 + (Sub10_diff_harvox - meanMatrix_harvox ) .^2 + (Sub11_diff_harvox - meanMatrix_harvox ) .^2 + (Sub12_diff_harvox - meanMatrix_harvox ) .^2 + (Sub13_diff_harvox - meanMatrix_harvox ) .^2 + (Sub15_diff_harvox - meanMatrix_harvox ) .^2 + (Sub17_diff_harvox - meanMatrix_harvox ) .^2 + (Sub18_diff_harvox - meanMatrix_harvox ) .^2 + (Sub19_diff_harvox - meanMatrix_harvox ) .^2 + (Sub20_diff_harvox- meanMatrix_harvox ) .^2;
std = sqrt(temp./15);
%std(abs(significant_harvox)<=0.007) = 1e-9;
std(abs(std)<=0.0126) = 1e-9;
colormap('parula')
imagesc(std)
colorbar



% flag the high average
significant_harvox = meanMatrix_harvox;
significant_harvox(abs(significant_harvox)<=0.007) = 1e-9;
colormap('parula')
imagesc(significant_harvox)
colorbar
title('Average Difference from Harvard-Oxford','FontSize',14)
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
temp = zeros(1,15);
for i = 1:48
    for j = 1:48
        lst = [Sub01_diff_harvox(i,j) Sub02_diff_harvox(i,j) Sub03_diff_harvox(i,j) Sub04_diff_harvox(i,j) Sub06_diff_harvox(i,j) Sub09_diff_harvox(i,j) Sub10_diff_harvox(i,j) Sub11_diff_harvox(i,j) Sub12_diff_harvox(i,j) Sub13_diff_harvox(i,j) Sub15_diff_harvox(i,j) Sub17_diff_harvox(i,j) Sub18_diff_harvox(i,j) Sub19_diff_harvox(i,j) Sub20_diff_harvox(i,j)];
        temp = cat(1,temp,lst);
    end
end
final = temp(2:end,:);

% find the indices of 20 flagged regions
linearIndexes = find(abs(significant_harvox) > 1e-9); 
[rows, columns] = ind2sub(size(significant_harvox), linearIndexes);


% 20 vs all
h_lst = zeros(48);
for i = 1:20
    r_ind = rows(i);
    temp = zeros(1,48*48);
    for j = 1:48*48
        [h,p] = ttest2(final((r_ind-1)*48 + columns(i),:), final(j,:));
        temp(j) = h;
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






