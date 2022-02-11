load_RPPA_data

%% How widespread is each drug effect? PCA
figure(2), clf
myDrug = "Dorsomorphin"; 
myDrug = "Ipatasertib";
myDrugs = {"AICAR";"CHIR-98014";"Cisplatin";"DMF";"DMSO";"DMSO*";"Dorsomorphin";"Everolimus";"Ipatasertib";"PF-00562271";"PF-4708671";"QNZ";"Stattic";"TAK-632";"U0126"};
myDrugs = {"AICAR";"CHIR-98014";"Cisplatin";"Dorsomorphin";"Everolimus";"Ipatasertib";"PF-00562271";"PF-4708671";"QNZ";"Stattic";"TAK-632";"U0126"};
%
for i=1:length(myDrugs)
figure(2), clf
myDrug = myDrugs{i};
%mySelection0 = tableCLwIC50(tableCLwIC50.Treatment==myDrug,:)
mySelection0 = dataFoldChangewIC50(dataFoldChangewIC50.Treatment==myDrug,:)
%mySelection0(mySelection0.CellLine=="MCF10A",:)=[];
mySelection0{:,3:end-1} = log2(mySelection0{:,3:end-1});
%mySelection0.AMPK_T172_(mySelection0.AMPK_T172_<-3)=-3;
mySelection0.AKT_T308_;
xaxis = 'IC50';

%Impute INFs
tmp = (mySelection0{:,3:(end-1)});
tmp(isinf(tmp)) = max(max(tmp(~isinf(tmp))));
tmp(isnan(tmp)) = 1;
mySelection0{:,3:(end-1)} = tmp;

%limit at 3
%tmp = (mySelection0{:,3:(end-1)});
%tmp(tmp>3) = 3;
%mySelection0{:,3:(end-1)} = tmp;

%myAnalytes  = mySelection0.Properties.VariableNames(3:(end-1))
myAnalytes = {'AKT_S473_' 'AKT_T308_' 'AMPK_T172_' 'cRAF_S338_' 'GSK3B_S9_' 'MAPK_T202_' 'MEK1_2_S217_' 'MET_T1234_' 'PDK1_S241_' 'SRC_Y527_' 'p38MAPK_T180_' 'gab1_Y627_' 'FAK_Y925_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'p53' 'CHK1_S345_' 'p27_T157_' 'PKCalpha_S657_' 'rsk_s380_' 'stat3_y750_' 'egfr_y992_' 'mtor_s2481_' 'egfr_y1068_' 'src_fam_y416_' 'bcl2_s70_' 'her2_y1248_' 'bcl_t56_' 'P70_S6K_T389_'}
%myAnalytes = {'p38MAPK_T180_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'CHK1_S345_' 'mtor_s2481_' 'P70_S6K_T389_'}
mask = ismember(mySelection0.Properties.VariableNames, [myAnalytes 'CellLine' 'Treatment', 'IC50']);
mySelection0(:,~mask) = [];

%PCA
%[wcoeff,score,latent,tsquared,explained] = pca(mySelection0{:,3:(end-1)}')
%scatter(score(:,1),score(:,2))%scatter(score(:,1),score(:,2), [], -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), 'Filled')
[wcoeff,score,latent,tsquared,explained] = pca(mySelection0{:,3:(end-1)})
scatter(score(:,1),score(:,2), [], -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), 'Filled')
CM = brewermap(3,'*Set1')
colormap(CM)
xlabel('Component 1')
ylabel('Component 2')
title(myDrug)
box on

figure(3)
%biplot(wcoeff(:,1:2),'Scores',score(:,1:2), 'VarLabels', mySelection0.CellLine)
biplot(wcoeff(:,1:2),'Scores',score(:,1:2),'Varlabels',mySelection0.Properties.VariableNames(3:(end-1)) );
pause
end

%% PCA of CHIR could disntinguish between sensitive and resistant - use only significant fold-changes

figure(12), clf
myDrug = "CHIR-98014";
%use only significant changes
mySelection0 = dataFoldChangewIC50(dataFoldChangewIC50.Treatment==myDrug,:)
%use all values
%mySelection0 = dataFoldChange_allwIC50(dataFoldChange_allwIC50.Treatment==myDrug,:)
%mySelection0(mySelection0.CellLine=="MCF10A",:)=[];
mySelection0{:,3:end-1} = log2(mySelection0{:,3:end-1});
xaxis = 'IC50';

%Impute INFs
tmp = (mySelection0{:,3:(end-1)});
tmp(isinf(tmp)) = max(max(tmp(~isinf(tmp))));
tmp(isnan(tmp)) = 1;
mySelection0{:,3:(end-1)} = tmp;

% %crop at |x|>3
% tmp = (mySelection0{:,3:(end-1)});
% limit = 3;
% tmp(tmp>limit) = limit;
% tmp(tmp<-limit) = -limit;
% mySelection0{:,3:(end-1)} = tmp;


pattern1 = {'p38MAPK_T180_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'CHK1_S345_' 'mtor_s2481_' 'P70_S6K_T389_'};
myAnalytes = {'AKT_S473_' 'AKT_T308_' 'AMPK_T172_' 'cRAF_S338_' 'GSK3B_S9_' 'MAPK_T202_' 'MEK1_2_S217_' 'MET_T1234_' 'PDK1_S241_' 'SRC_Y527_' 'p38MAPK_T180_' 'gab1_Y627_' 'FAK_Y925_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'p53' 'CHK1_S345_' 'p27_T157_' 'PKCalpha_S657_' 'rsk_s380_' 'stat3_y750_' 'egfr_y992_' 'mtor_s2481_' 'egfr_y1068_' 'src_fam_y416_' 'bcl2_s70_' 'her2_y1248_' 'bcl_t56_' 'P70_S6K_T389_'}
%myAnalytes = pattern1;
mask = ismember(mySelection0.Properties.VariableNames, [myAnalytes 'CellLine' 'Treatment', 'IC50']);
mySelection0(:,~mask) = [];

%PCA
%[wcoeff,score,latent,tsquared,explained] = pca(mySelection0{:,3:(end-1)}')
%scatter(score(:,1),score(:,2))%scatter(score(:,1),score(:,2), [], -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), 'Filled')
[wcoeff,score,latent,tsquared,explained] = pca(mySelection0{:,3:(end-1)}, 'NumComponents', 3)
%scatter(score(:,1),score(:,2), [], -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), 'Filled')
ph = gscatter(score(:,1),score(:,2), -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), brewermap(3, '*Set1'))
%legend(ph, 'Senistive', 'MCF10A', 'Resistant', 'Location', 'northeastoutside')
CM = brewermap(3,'*Set1')
colormap(CM)
xlabel('Component 1')
ylabel('Component 2')
title([convertStringsToChars(myDrug) ', only sign. changes'])
box on

%Annotation:
C1 = wcoeff(:,1); C1 = C1/sum(abs(C1));
C2 = wcoeff(:,2); C2 = C2/sum(abs(C2));
[~, sortmask1] = sort(abs(C1), 'descend');
C1s = C1(sortmask1)
C1sLabel = myAnalytes(sortmask1);
%[~, sortmask2] = sort(abs(C2), 'descend');
[~, sortmask2] = sort(C2, 'descend');
C2s = C2(sortmask2)
C2sLabel = myAnalytes(sortmask2);
hold on
xLim = xlim; dx = xLim(2)-xLim(1);
yLim = ylim; dy = yLim(2)-yLim(1);
plot([-0.0*dx 26], [0 0],'k--'), plot(26+0.025*dx*[-1 0 -1], 0.015*dy*[1 0 -1], 'k-')
text(18,0, C1sLabel(1:5), 'Interpreter', 'none', 'VerticalAlignment', 'top')
plot([0 0], [-0.0*dy 1.8],'k--'), plot(0.015*dx*[1 0 -1],1.8+0.025*dy*[-1 0 -1],'k-')
text(0.+0.02*dx,1.8, C2sLabel(1:5), 'Interpreter', 'none', 'VerticalAlignment', 'top')

legend(ph, 'Sensitive', 'MCF10A', 'Resistant', 'Location', 'northeastoutside')
fpos = get(gcf, 'Position');
fsize = [500 360];
set(gcf, 'Position', [fpos(1:2) fsize])
print(gcf, ['../figures//PCA_scores_' convertStringsToChars(myDrug) '_onlysign.png'], '-dpng', '-r600')

figure(13)
%biplot(wcoeff(:,1:2),'Scores',score(:,1:2), 'VarLabels', mySelection0.CellLine)
ph = biplot(wcoeff(:,1:2),'Scores',score(:,1:2),'Varlabels',mySelection0.Properties.VariableNames(3:(end-1)) );
box on
fpos = get(gcf, 'Position');
set(gcf, 'Position', [fpos(1:2) fsize])

print(gcf, ['../figures/PCA_biplot_' convertStringsToChars(myDrug) '_onlysign.png'], '-dpng', '-r600')


%% CHIR could disntinguish between sensitive and resistant - USE ALL
figure(2), clf
myDrug = "CHIR-98014";
%use only significant changes
%mySelection0 = dataFoldChangewIC50(dataFoldChangewIC50.Treatment==myDrug,:)
%use all values
mySelection0 = dataFoldChange_allwIC50(dataFoldChange_allwIC50.Treatment==myDrug,:)
%mySelection0(mySelection0.CellLine=="MCF10A",:)=[];
mySelection0{:,3:end-1} = log2(mySelection0{:,3:end-1});
xaxis = 'IC50';

%Impute INFs
tmp = (mySelection0{:,3:(end-1)});
tmp(isinf(tmp)) = max(max(tmp(~isinf(tmp))));
tmp(isnan(tmp)) = 1;
mySelection0{:,3:(end-1)} = tmp;

% %crop at |x|>3
% tmp = (mySelection0{:,3:(end-1)});
% limit = 3;
% tmp(tmp>limit) = limit;
% tmp(tmp<-limit) = -limit;
% mySelection0{:,3:(end-1)} = tmp;


pattern1 = {'p38MAPK_T180_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'CHK1_S345_' 'mtor_s2481_' 'P70_S6K_T389_'};
myAnalytes = {'AKT_S473_' 'AKT_T308_' 'AMPK_T172_' 'cRAF_S338_' 'GSK3B_S9_' 'MAPK_T202_' 'MEK1_2_S217_' 'MET_T1234_' 'PDK1_S241_' 'SRC_Y527_' 'p38MAPK_T180_' 'gab1_Y627_' 'FAK_Y925_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'p53' 'CHK1_S345_' 'p27_T157_' 'PKCalpha_S657_' 'rsk_s380_' 'stat3_y750_' 'egfr_y992_' 'mtor_s2481_' 'egfr_y1068_' 'src_fam_y416_' 'bcl2_s70_' 'her2_y1248_' 'bcl_t56_' 'P70_S6K_T389_'}
%myAnalytes = pattern1;
mask = ismember(mySelection0.Properties.VariableNames, [myAnalytes 'CellLine' 'Treatment', 'IC50']);
mySelection0(:,~mask) = [];

%PCA
%[wcoeff,score,latent,tsquared,explained] = pca(mySelection0{:,3:(end-1)}')
%scatter(score(:,1),score(:,2))%scatter(score(:,1),score(:,2), [], -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), 'Filled')
[wcoeff,score,latent,tsquared,explained] = pca(mySelection0{:,3:(end-1)}, 'NumComponents', 3)
%scatter(score(:,1),score(:,2), [], -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), 'Filled')
ph = gscatter(score(:,1),score(:,2), -1*(mySelection0.IC50<2)+(mySelection0.IC50>6), brewermap(3, '*Set1'))
%legend(ph, 'Senistive', 'MCF10A', 'Resistant', 'Location', 'northeastoutside')
CM = brewermap(3,'*Set1')
colormap(CM)
xlabel('Component 1')
ylabel('Component 2')
title(myDrug)
box on

%Annotation:
C1 = wcoeff(:,1); C1 = C1/sum(abs(C1));
C2 = wcoeff(:,2); C2 = C2/sum(abs(C2));
[~, sortmask1] = sort(abs(C1), 'descend');
C1s = C1(sortmask1)
C1sLabel = myAnalytes(sortmask1);
%[~, sortmask2] = sort(abs(C2), 'descend');
[~, sortmask2] = sort(C2, 'descend');
C2s = C2(sortmask2)
C2sLabel = myAnalytes(sortmask2);
hold on
yLim = ylim; ylim([yLim(1) 1.1*yLim(2)]);
xLim = xlim; dx = xLim(2)-xLim(1);
yLim = ylim; dy = yLim(2)-yLim(1);
plot([-0.0*dx 26], [0 0],'k--'), plot(26+0.025*dx*[-1 0 -1], 0.015*dy*[1 0 -1], 'k-')
text(18,0, C1sLabel(1:5), 'Interpreter', 'none', 'VerticalAlignment', 'top')
plot([0 0], [-0.0*dy 2.8],'k--'), plot(0.015*dx*[1 0 -1],2.8+0.025*dy*[-1 0 -1],'k-')
text(0.+0.02*dx,2.8, C2sLabel(1:6), 'Interpreter', 'none', 'VerticalAlignment', 'top')

legend(ph, 'Sensitive', 'MCF10A', 'Resistant', 'Location', 'northeastoutside')
fpos = get(gcf, 'Position');
fsize = [500 360];
set(gcf, 'Position', [fpos(1:2) fsize])

print(gcf, ['../figures/PCA_scores_' convertStringsToChars(myDrug) '.png'], '-dpng', '-r600')

figure(3)
%biplot(wcoeff(:,1:2),'Scores',score(:,1:2), 'VarLabels', mySelection0.CellLine)
ph = biplot(wcoeff(:,1:2),'Scores',score(:,1:2),'Varlabels',mySelection0.Properties.VariableNames(3:(end-1)) );
box on
fpos = get(gcf, 'Position');
set(gcf, 'Position', [fpos(1:2) fsize])

print(gcf, ['../figures/PCA_biplot_' convertStringsToChars(myDrug) '.png'], '-dpng', '-r600')




