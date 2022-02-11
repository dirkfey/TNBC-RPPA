%% Back it up by plotting the Pattern mTOR p70S6K CHK1 p38 S6RIB in sensitive vs resistant
myAnalytes = {'p38MAPK_T180_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'CHK1_S345_' 'mtor_s2481_' 'P70_S6K_T389_'}
myDrugs = ["AICAR";"CHIR-98014";"Dorsomorphin";"QNZ";"TAK-632"];
nDrugs = length(myDrugs);
figure(24000), clf
subplot1(1,nDrugs+1)
for i=1:nDrugs
    myDrug = myDrugs(i)
    %use this to omit non-significant changes (0-impute)
    %mySelection0 = dataFoldChangewIC50(dataFoldChangewIC50.Treatment==myDrug,:);
    %use this to use original data
    mySelection0 = dataFoldChange_allwIC50(dataFoldChange_allwIC50.Treatment==myDrug,:);
    yy_sign = dataFCSignificant(dataFoldChange_allwIC50.Treatment==myDrug,:);
    %
    [mySelection00, idx] = sortrows(mySelection0, 'IC50');
    yy_sign = yy_sign(idx,:);
    figure(23000+i), clf
    %plot(log10(mySelection00.IC50), mySelection00{:,convertStringsToChars(myAnalytes)},'o:')
    xx = 1:size(mySelection00,1);
    xx = [0.8 1 1.2 2 2.8 3 3.2];
    yy = mySelection00{:,convertStringsToChars(myAnalytes)};
    yy_sign = yy_sign{:,convertStringsToChars(myAnalytes)}; 
    yy = log2(yy);
    ph = plot(xx, yy,'o:')
    hold on
    %mark significant changes:
    set(gca, 'ColorOrderIndex', 1)
    %set(gcf, 'ColorOrderIndex', 1)
    removematrix = yy_sign; %removematrix(yy_sign==0)= -1;
    xxx = repmat(xx', 1, size(yy,2));
    plot(xxx.*(removematrix), yy,'.', 'MarkerSize', 18);
    set(gca, 'xLim', [0.5 3.5])
    %
    set(gca, 'xTick', [1 2 3], 'XTickLabel', {'Sens.', 'MCF10A', 'Resist.'})
    %plot(log10(mySelection00.IC50), mean(mySelection00{:,convertStringsToChars(myAnalytes)},2),'ko-')
    plot(xx, mean(yy,2),'ko-')
    %change FC1 to NAN for mean calculation
    mySelection000=yy;
    mySelection000(mySelection000==1)=nan;
    %plot(xx, nanmean(mySelection000,2),'rx--')
    %
    %Stats:
    %test each analyte seperately
    [h, p] = ttest2(yy(1:3,:),yy(5:7,:));
    %test cluster: lump and test all analytes together
    [hm, pm] = ttest2(reshape(yy(1:3,:),[],1),reshape(yy(5:7,:),[],1));
    h_plotpval(pm)
    %Do also Mann-Whitney U-test (Wilcoxon)
    %[pm, hm] = ranksum(reshape(yy(1:3,:),[],1),reshape(yy(5:7,:),[],1))
    %h_plotpval(pm)
    %highlight most significant individual analyte
    [minp idx] = min(p);
    set(ph(idx), 'LineStyle', '-')
    %
    title(myDrug)
    legend(ph, myAnalytes, 'Location','eastoutside')
    %Save figure:
    %print(gcf, ['./figs/pattern1_all_' char(myDrug) '.png'], '-dpng', '-r600')
    %
    % Aggregated DISTRIBUTIONPLOTS using original data (pval<0.5 value instead of 0)
    figure(24000)
    mySelection0 = dataFoldChange_allwIC50(dataFoldChange_allwIC50.Treatment==myDrug,:);
    mySelection00 = sortrows(mySelection0, 'IC50')
    %old: notsign = yy==0;
    notsign = ~yy_sign;
    yy = mySelection00{:,convertStringsToChars(myAnalytes)};
    yy = log2(yy);
    %subplot1(1,6)
    subplot1(i)
    catcolors = brewermap(7,'dark2');
    tmp = [reshape(yy(1:3,:),[],1),reshape(yy(5:7,:),[],1)];
    notsign = [reshape(notsign(1:3,:),[],1),reshape(notsign(5:7,:),[],1)];
    psh1 = plotSpread(tmp, 'categoryIdx' , [reshape(repmat(myAnalytes,3,1),[],1) reshape(repmat(myAnalytes,3,1),[],1)], 'categoryColors', catcolors);
    %plotSpread(tmp, 'categoryIdx' , [reshape(repmat([1:7],3,1),[],1) reshape(repmat([1:7],3,1),[],1)]);
    psh2 = plotSpread(tmp, 'categoryIdx' , notsign, 'categorymarkers', {'o' 'o'});
    %set(psh{1}(2,1), 'MarkerSize', 15)
    set(psh1{1}(:,:), 'MarkerSize', 16)
    set(psh2{1}(:,1), 'MarkerEdgeColor', 'none')
    set(psh2{1}(:,2), 'MarkerEdgeColor', [1 1 1], 'Linewidth', 1, 'MarkerSize', 4)
    %annotate pval:
    [h, pval] = ttest2(tmp(:,1), tmp(:,2));
    if pval<0.001
        pval = sprintf('p-val=%1.1e', pval);
    else
        pval = sprintf('p-val=%1.3f', pval);
    end
    text(1.6, 1.5, pval, 'HorizontalAlignment', 'center')
    %inlcude Boxplot:
    hold on
    boxplot(tmp, 'notch', 'on', 'Symbol', '', 'Whisker', 0, 'Colors', 'k' )
    set(findobj(gca, 'Tag', 'Median'), 'LineWidth', 2)
    set(gca, 'ylim', [-2 2])
    set(gca, 'xTick', [1 2], 'XTickLabel', {'Sens.' 'Resist.'})
    title(myDrug)
end
%
fpos = get(gcf, 'Pos');
set(gcf, 'Position', [fpos(1:2) 1200 260])

lh = legend(psh1{1}(1,:), myAnalytes, 'Interpreter', 'none', 'Location', 'northeast');
lpos = get( lh, 'Pos');
set(lh, 'Pos', [0.81 lpos(2:4)])
fh = get(gcf, 'Children');
delete(fh(1))


axes('Position', [0.81 0.2 lpos(3) 0.15], 'xTick', [], 'yTick', [], 'YDir', 'reverse', 'YLim', [0 2])

box on
plot(0.15, 0.5, '.', 'Markersize', 16, 'Color', 0.33*[1 1 1])
text(0.3, 0.5, 'p-val<0.05')
hold on
plot(0.15, 1.5, '.', 'Markersize', 8, 'Color', 0.33*[1 1 1])
text(0.3, 1.5, 'not significant')
text(0.01, 2.5, '4 replicates, two-sample')
text(0.01, 3.3, 't-test')
title('Fold change significance')
set(gca,'yLim', [0 2], 'xlim', [0, 1], 'xTick', [], 'yTick', [], 'YDir', 'reverse')

print(gcf, ['../figures/boxplot_pattern1_tmp.png'], '-dpng', '-r600')
