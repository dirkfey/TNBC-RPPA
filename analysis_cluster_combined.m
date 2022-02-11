

%% Clustergramms for all drugs - Heatmaps

myDrugs = ["AICAR";"CHIR-98014";"Cisplatin";"Dorsomorphin";"Everolimus";"Ipatasertib";"PF-00562271";"PF-4708671";"QNZ";"Stattic";"TAK-632";"U0126"];
%myDrugs = ["AICAR";"CHIR-98014";"Ipatasertib"];

TF = contains(dataFoldChangewIC50.Treatment, myDrugs);
mySelection0 = dataFoldChange_allwIC50(TF,:);
mySelection0 = dataFoldChangewIC50(TF,:);
mySelection0{:,3:end-1} = log2(mySelection0{:,3:end-1});
xaxis = 'IC50';

mySelection0 = sortrows(mySelection0,"IC50");
mySelection0 = sortrows(mySelection0,"Treatment")

%Impute INFs
tmp = (mySelection0{:,3:(end-1)});
tmp(isinf(tmp)) = max(max(tmp(~isinf(tmp))));
tmp(isnan(tmp)) = 1;
mySelection0{:,3:(end-1)} = tmp;

myAnalytes = {'AKT_S473_' 'AKT_T308_' 'AMPK_T172_' 'cRAF_S338_' 'GSK3B_S9_' 'MAPK_T202_' 'MEK1_2_S217_' 'MET_T1234_' 'PDK1_S241_' 'SRC_Y527_' 'p38MAPK_T180_' 'gab1_Y627_' 'FAK_Y925_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'p53' 'CHK1_S345_' 'p27_T157_' 'PKCalpha_S657_' 'rsk_s380_' 'stat3_y750_' 'egfr_y992_' 'mtor_s2481_' 'egfr_y1068_' 'src_fam_y416_' 'bcl2_s70_' 'her2_y1248_' 'bcl_t56_' 'P70_S6K_T389_', ...
    'PTEN', 'PI3KalPha'};
%myAnalytes = {'AKT_S473_' 'AKT_T308_' 'PDK1_S241_' 'p27_T157_' , ...
%    'PTEN', 'PI3KalPha'};
mask = ismember(mySelection0.Properties.VariableNames, [myAnalytes 'CellLine' 'Treatment', 'IC50']);
mySelection0(:,~mask) = [];

%%
%Convert to DataMatrix and cluster
import bioma.data.*
mySelection00 = mySelection0;
rownames = join([mySelection00.CellLine, mySelection00.Treatment],'::')
%rownames = join([mySelection00.Treatment, mySelection00.CellLine],'::')
colordefine = cell(size(rownames));
colordefine(mySelection00.IC50<5)={[0.5 0.5 1]};
colordefine(mySelection00.IC50>5)={[0.8 0.8 0.8]};
colordefine(mySelection00.IC50>6)={[1 0.5 0.5]};
s = struct('Labels', rownames, 'Colors', colordefine);
myDM = DataMatrix(mySelection00{:,3:end},'RowNames', rownames,'ColNames', mySelection00.Properties.VariableNames(3:end) )
cgo = clustergram(myDM(:,1:(end-1)), 'cluster', 2, 'Colormap', 'redbluecmap')
%set(cgo, 'RowLabels', rownames, 'LabelsWithMarkers', 1, 'RowLabelsColor', s)
th = addTitle(cgo, myDrug)
%set(cgo, 'Fontsize', 10)
figure(3), clf
fh = plot(cgo,gcf);
%ha = get(fh, 'Children'); ha = ha(2);
%set(fh, 'Position', [apos(1) 0.3 apos(3), 0.7*apos(4)])
%fh.Title.String = myDrug;
%tpos = 8;
%fh.Title.Position(2)=1.4*tpos;
%print(gcf, ['../figures/cluster_2_' char(myDrug) '.png'], '-dpng', '-r600')

apos = get(fh, 'Position')

%%
i=1;
for i = 1%1:length(myDrugs)
figure(3), clf
th = addTitle(cgo, myDrugs(i))
th.String = myDrugs(i)
fh = plot(cgo,gcf);
box on
fpos = get(gcf, 'Pos');
set(gcf, 'Pos', [fpos(1:2) 560 300])
set(fh, 'ylim', (i-1)*7+[0.5 7.5])
set(fh, 'Position', [apos(1) 0.37 apos(3) apos(4)+apos(2)-0.37])
set(fh, 'Fontsize', 10)
drawnow
print(gcf, ['../figures/cluster_2_combined_' char(myDrugs(i)) '.png'], '-dpng', '-r600')
end