%% Analysis of responsiveness:
%  Summarize number of significant changes (in netwok nodes) per cell line
%  (summarised over treatments: number; average over treatments: fold-change (FC))

%% Select drugs and analytes and pre-process data:

myDrugs = ["AICAR";"CHIR-98014";"Cisplatin";"Dorsomorphin";"Everolimus";"Ipatasertib";"PF-00562271";"PF-4708671";"QNZ";"Stattic";"TAK-632";"U0126"];
myAnalytes = {'AKT_S473_' 'AKT_T308_' 'AMPK_T172_' 'cRAF_S338_' 'GSK3B_S9_' 'MAPK_T202_' 'MEK1_2_S217_' 'MET_T1234_' 'PDK1_S241_' 'SRC_Y527_' 'p38MAPK_T180_' 'gab1_Y627_' 'FAK_Y925_' 'mTOR_S2448_' 'S6RIB_S240_' 'S6RIB_S235_' 'p53' 'CHK1_S345_' 'p27_T157_' 'PKCalpha_S657_' 'rsk_s380_' 'stat3_y750_' 'egfr_y992_' 'mtor_s2481_' 'egfr_y1068_' 'src_fam_y416_' 'bcl2_s70_' 'her2_y1248_' 'bcl_t56_' 'P70_S6K_T389_', ...
    'PTEN', 'PI3KalPha'};

TF = contains(dataFCSignificant.treatments, myDrugs);
dataFCSign_selected = dataFCSignificant(TF,[{'cellLines' 'treatments'} myAnalytes])
%
dataFC_selected = dataFoldChange(TF,[{'cellLines' 'treatments'} myAnalytes]);
% replace NaN (here zeros) with 1:
dataFC_selected{:,myAnalytes}(dataFC_selected{:,myAnalytes}==0)=1; 
% log2 transform:
dataFC_selected{:, myAnalytes} = log2(dataFC_selected{:,myAnalytes})

%% Compute number of responseive network nodes per treatment and summarize result
dataFCSign_selected.Number_responsive = sum(dataFCSign_selected{:,myAnalytes},2);
responsive = sortrows ( ...
  groupsummary(dataFCSign_selected,'cellLines','sum','Number_responsive'), ...
  'sum_Number_responsive', 'descend')
%% Compute the mean average FC per analyte and treatment 
dataFC_selected.AvgFC = mean(dataFC_selected{:,myAnalytes},2)
responsiveFC = sortrows ( ...
  groupsummary(dataFC_selected,'cellLines','mean','AvgFC'), ...
  'mean_AvgFC', 'descend')

