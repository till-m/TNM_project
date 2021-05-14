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



%% Analysis
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
caxis([-0.05 0.05]) 
%set(gca,'xtick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)

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
set(gca,'xtick',[1:48])
set(gca,'ytick',[1:48])
%set(gca,'ytick',[1 round(size(output.Ep.A,1)/2) size(output.Ep.A,1)])
xlabel('region (from)','FontSize',12)
ylabel('region (to)','FontSize',12)
