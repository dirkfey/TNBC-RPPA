

%% Clustergramms for all drugs - Heatmaps

myDrugs = {"AICAR";"CHIR-98014";"Cisplatin";"Dorsomorphin";"Everolimus";"Ipatasertib";"PF-00562271";"PF-4708671";"QNZ";"Stattic";"TAK-632";"U0126"};
myDrugs = {"Ipatasertib";};

for i=1:length(myDrugs)
myDrug = myDrugs{i};
%mySelection0 = tableCLwIC50(tableCLwIC50.Treatment==myDrug,:)
mySelection0 = dataFoldChangewIC50(dataFoldChangewIC50.Treatment==myDrug,:)
mySelection0 = dataFoldChange_allwIC50(dataFoldChangewIC50.Treatment==myDrug,:)
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

myAnalytes = {'AKT_S473_' 'AKT_T308_' 'AMPK_T172_' 'cRAF_S338_' 'GSK3B_S9_' 'MAPK_T202_' 'MEK1_2_S217_' 'MET_T1234_' 'PDK1_S241_' 'SRC_Y527_' 'p38MAPK_T180_' 'gab1_Y627_' 'FAK_Y925_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'p53' 'CHK1_S345_' 'p27_T157_' 'PKCalpha_S657_' 'rsk_s380_' 'stat3_y750_' 'egfr_y992_' 'mtor_s2481_' 'egfr_y1068_' 'src_fam_y416_' 'bcl2_s70_' 'her2_y1248_' 'bcl_t56_' 'P70_S6K_T389_', ...
    'PTEN', 'PI3KalPha'};
myAnalytes = {'AKT_S473_' 'AKT_T308_' 'PDK1_S241_' 'p27_T157_' , ...
    'PTEN', 'PI3KalPha'};
mask = ismember(mySelection0.Properties.VariableNames, [myAnalytes 'CellLine' 'Treatment', 'IC50']);
mySelection0(:,~mask) = [];

%Convert to DataMatrix and cluster
import bioma.data.*
mySelection00 = sortrows(mySelection0, 'IC50')
myDM = DataMatrix(mySelection00{:,3:end},'RowNames', mySelection00.CellLine,'ColNames', mySelection00.Properties.VariableNames(3:end) )
cgo = clustergram(myDM(:,1:(end-1)), 'cluster', 2, 'Colormap', 'redbluecmap')
%addTitle(cgo, myDrug)
fh = plot(cgo);
%ha = get(fh, 'Children'); ha = ha(2);
apos = get(fh, 'Position');
%set(fh, 'Position', [apos(1) 0.3 apos(3), 0.7*apos(4)])
fh.Title.String = myDrug;
tpos = 8;
fh.Title.Position(2)=1.4*tpos;
%print(gcf, ['../figures/cluster_2_' char(myDrug) '.png'], '-dpng', '-r600')
end




