%% Differential expression analysis of (phospho-)proteins: resistant vs sensitive cells
%  Requires: dataMean

%load dataRPPAmean.mat
load mat_dataMean.mat

%% Pre-process data:

data_ctr = dataMean(dataMean.Treatment == "DMSO",:);

%Q: What about DMSO*?
data_ctr = dataMean(ismember(dataMean.Treatment,["DMSO" "DMF"]),:);

resistantcells = ["MFM223" "MDA-MB-157" "HDQP-1"];
sensitivecells = ["MDA-MB-468" "CAL-85-1" "MDA-MB-231"];

data_res = data_ctr(ismember(data_ctr.CellLine, resistantcells),:)
data_sen = data_ctr(ismember(data_ctr.CellLine, sensitivecells),:)

import bioma.data.*
expr_res = DataMatrix(data_res{:,3:end}', 'Colnames', data_res.CellLine, 'RowNames', data_res.Properties.VariableNames(3:end));
expr_sen = DataMatrix(data_sen{:,3:end}', 'Colnames', data_sen.CellLine, 'RowNames', data_sen.Properties.VariableNames(3:end));
%log transform:
expr_res = log2(expr_res);
expr_sen = log2(expr_sen);

%% Perform differential analysis using permutation two-sample t-tests
[pvalues,tscores] = mattest(expr_sen, expr_res, ...
                    'permute',1000);
[pFDR, qvalues] = mafdr(pvalues);

clear diffstruct upstruct downstruct
diffstruct = mavolcanoplot(expr_sen, expr_res, pFDR, 'PCutoff', 0.05, 'FoldChange', 2)
diffstruct_fc1p5 = mavolcanoplot(expr_sen, expr_res, pFDR, 'PCutoff', 0.05, 'FoldChange', 1.5)

%% Visualise results #1: Plotspead and Boxplot figures of |log2 FC| > 2 & pval < 0.05
for i=1:length(diffstruct.PValues)
    analyte = diffstruct.PValues.RowNames{i};
    pval = diffstruct.PValues(analyte,1);
    fc = diffstruct.FoldChanges(analyte,1);
    figure(i), clf
    tmp = [expr_sen(analyte,:) expr_res(analyte,:)]
    %plotSpread(double([expr_sen(analyte,:); expr_res(analyte,:)]'))
    mydata = double([expr_sen(analyte,:) expr_res(analyte,:)]);
    myG = [repmat({'sen'}, 1, 6) repmat({'res'}, 1, 6)];
    bph  = boxplot(mydata', myG', 'Widths', 0.9, 'Whisker', 0 );
    psh = plotSpread(mydata', 'distributionIdx', myG')
    %
    set(psh{1}, 'MarkerSize', 10, 'Color', [0 0 0.8])
    set(bph(7,:), 'Marker', 'none' )   
    set(bph(5,:), 'Color', 0.5*[1 1 1])   
    set(bph(6,:), 'Color', [0.8 0 0], 'LineWidth', 1)   
    %
    ylabel([analyte ' (log2)'], 'Interpreter', 'none' )
    title(sprintf('FC = %1.1f, p = %1.1e', fc, pval), 'Interpreter', 'none')
    fpos = get(gcf, 'Position');
    set(gcf, 'Position', [fpos(1:2), 220 200])
    %
    print(sprintf('../figures/fig_sen_vs_res_boxplot_%u_%s.png', i, analyte), '-dpng', '-r600')
end

%% Visualise results #1: Look at PI3K/AKT module in particular: |log2 FC| > 2 & pval < 0.05

% Focus on these analytes:
myanalytes = {'PI3KalPha', 'PDK1', 'PDK1_S241_', 'AKT_S473_', 'p27_T157_'};
% Make distribution boxplots:
for i=1:length(myanalytes)
    analyte = myanalytes{i};
    pval = diffstruct_fc1p5.PValues(analyte,1);
    fc = diffstruct_fc1p5.FoldChanges(analyte,1);
    figure(i), clf
    tmp = [expr_sen(analyte,:) expr_res(analyte,:)];
    %plotSpread(double([expr_sen(analyte,:); expr_res(analyte,:)]'))
    mydata = double([expr_sen(analyte,:) expr_res(analyte,:)]);
    myG = [repmat({'sen'}, 1, 6) repmat({'res'}, 1, 6)];
    bph  = boxplot(mydata', myG', 'Widths', 0.9, 'Whisker', 0 );
    psh = plotSpread(mydata', 'distributionIdx', myG');
    %
    set(psh{1}, 'MarkerSize', 10, 'Color', [0 0 0.8])
    set(bph(7,:), 'Marker', 'none' )   
    set(bph(5,:), 'Color', 0.5*[1 1 1])   
    set(bph(6,:), 'Color', [0.8 0 0], 'LineWidth', 1)   
    %
    ylabel([analyte ' (log2)'], 'Interpreter', 'none' )
    title(sprintf('FC = %1.1f, p = %1.1e', fc, pval), 'Interpreter', 'none')
    fpos = get(gcf, 'Position');
    set(gcf, 'Position', [fpos(1:2), 220 200])
    %
    print(sprintf('../figures/fig_sen_vs_res_boxplot_AKTaxis_%u_%s.png', i, analyte), '-dpng', '-r600')
end

%% Visualise results #4: Heatmap clustergrams 
expr = [expr_sen expr_res];
%expr('Casp7_cl',:) = [];
[imputeidx1, imputeidx2]= find(expr<1e-6);
expr(imputeidx1, imputeidx2)=NaN;
expr = expr-nanmean(expr,2);
expr(imputeidx1, imputeidx2)=0;
% Clustergramm #1 for |log2 FC| > 2 & pval < 0.05:
cgo1 = clustergram(expr(diffstruct.GeneLabels,:), 'Cluster', 1, 'Colormap', 'redbluecmap');
% Clustergramm #1 for |log2 FC| > 1.5 & pval < 0.05:
cgo2 = clustergram(expr(diffstruct_fc1p5.GeneLabels,:), 'Cluster', 1, 'Colormap', 'redbluecmap');

%% Export clustergrams to heatmap figure
plot(cgo1)
fpos = get(gcf, 'Position');
set(gcf,'Position', [fpos(1:2) 450 650])

print(gcf, '../figures/heatmap_sen_vs_res.svg', '-dsvg')
print(gcf, '../figures/heatmap_sen_vs_res.png', '-dpng', '-r600')

plot(cgo2)
fpos = get(gcf, 'Position');
set(gcf,'Position', [fpos(1:2) 450 850])

print(gcf, '../figures/heatmap_sen_vs_res_fc1p5.svg', '-dsvg')
print(gcf, '../figures/heatmap_sen_vs_res_fc1p5.png', '-dpng', '-r600')
