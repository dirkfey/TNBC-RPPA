
%load_RPPA_data ... skipt that imports the data from ".\rawdata\TNBC RPPA Panel 29.4.19 clean.xlsx"  
%   and processes it:
%   1. Filter, 
%   2. Calculate averages over replicates
%   3. Calculate mean-fold changes and test statistical significance 
%
%   Input: .\rawdata\TNBC RPPA Panel 29.4.19 clean.xlsx
%   Outputs: 
%   1. dataRPPAraw_.mat, here
%   2. mat_dataMean.mat, here
%      .\processeddata\RPPA_mean.xlsx
%   3. mat_dataFoldChange.mat 
%      ./processeddata/dataFoldChange.xlsx
%   Dependencies: none

dataRPPAraw = readtable('.\rawdata\TNBC RPPA Panel 29.4.19 clean.xlsx')

%Remove lines 94, 96, 112

dataRPPAraw([94, 96, 112]-1, :) = [];
%%
dataIC50 = table(...
    {'CAL-85-1';'HDQP-1';'MCF10A';'MDA-MB-157';'MDA-MB-231';'MDA-MB-468';'MFM223'},... 
    [0.4;       6.4;        5.1;    8.0;         1.9;         0.2;             10],...
    'Rownames', ...
    {'CAL_85_1';'HDQP_1';'MCF10A';'MDA_MB_157';'MDA_MB_231';'MDA_MB_468';'MFM223'},...
    'Variablenames', {'CellLine' 'IC50'});

%% remove datapoints 1e4 
dataRPPAfilt = dataRPPAraw(1:(end-2),:); 
data = dataRPPAfilt{:,4:end};
data(data>1e4) = nan; 
dataRPPAfilt{:,4:end} = data;

%%
save dataRPPAraw_.mat dataRPPAraw dataRPPAfilt dataIC50

%% Average over replicates: calculate mean for each cellline and treatment 
[G, cellline, treatment] = findgroups(dataRPPAfilt.CellLine,dataRPPAfilt.Treatment)

mymean = @(x) nanmean(x,1); 

meandata = splitapply(mymean,dataRPPAfilt{:,4:end},G)
data_DMSO = dataRPPAfilt(dataRPPAfilt.Treatment=="DMSO",:)

dataMean = [cellline, treatment, array2table( meandata)];
dataMean.Properties.VariableNames = dataRPPAfilt.Properties.VariableNames([1 3:end])

dataMeanwIC50 = join(dataMean, dataIC50, 'LeftKeys', 'CellLine','RightKeys', 'CellLine');

%%
writetable(dataMean, '.\processeddata\RPPA_mean.xlsx', 'Sheet', 'RPPA mean')
writetable(table('Source:','load_RPPA_data.m', cd), '.\processeddata\RPPA_mean.xlsx', 'Sheet', 'anno')
save mat_dataMean.mat dataMean

%% Test signficance of treatment changes
[G, cellLines, treatments] = findgroups(dataRPPAfilt.CellLine, dataRPPAfilt.Treatment);
varNames = dataRPPAfilt.Properties.VariableNames(4:end);
Nvars = size(dataRPPAfilt(1,4:end),2);
%init:
changefold = nan(max(G), Nvars);
changefold_all = changefold;
changepval = changefold;
changesign = changefold;
for i=1:max(G)
    i_dataTreatment = dataRPPAfilt{G==i,4:end};
    i_dataDMSO = dataRPPAfilt{strcmp(dataRPPAfilt.CellLine,cellLines{i})&(dataRPPAfilt.Treatment=="DMSO"),4:end};
    %
    %dataRPPAfilt(G==i,1:end)
    %dataRPPAfilt(strcmp(dataRPPAfilt.CellLine,cellLines{i})&(dataRPPAfilt.Treatment=="DMSO"),1:end)
    %
    %mytest = @(X) ttest(X,);
    %keyboard
    for ii = 1:Nvars
        [flag, pvalue] = ttest2(i_dataTreatment(:,ii),i_dataDMSO(:,ii));
        changesign(i,ii) = flag;
        changepval(i,ii) = pvalue;
        changefold_all(i,ii) = mean(i_dataTreatment(:,ii),'omitnan')./mean(i_dataDMSO(:,ii), 'omitnan');
        if isnan(flag), flag = 0; end
        if flag
            changefold(i,ii) = mean(i_dataTreatment(:,ii),'omitnan')./mean(i_dataDMSO(:,ii), 'omitnan');
        else
            changefold(i,ii) = 1;
        end
    end
end
%%
dataFoldChange = [table(cellLines, treatments) array2table(changefold, 'VariableNames', varNames)];
dataFoldChange_all = [table(cellLines, treatments) array2table(changefold_all, 'VariableNames', varNames)];
dataFCSignificant = [table(cellLines, treatments) array2table(changesign, 'VariableNames', varNames)];
dataFCpvalue = [table(cellLines, treatments) array2table(changepval, 'VariableNames', varNames)];

dataFoldChangewIC50 = join(dataFoldChange, dataIC50, 'LeftKeys', 'cellLines','RightKeys', 'CellLine');
dataFoldChangewIC50.Properties.VariableNames(1:2) = {'CellLine', 'Treatment'};
dataFoldChange_allwIC50 = join(dataFoldChange_all, dataIC50, 'LeftKeys', 'cellLines','RightKeys', 'CellLine');
dataFoldChange_allwIC50.Properties.VariableNames(1:2) = {'CellLine', 'Treatment'};
head(dataFoldChangewIC50)

%%
writetable(dataFoldChange_allwIC50, './processeddata/dataFoldChange.xlsx', 'Sheet', 'FC all')
writetable(dataFoldChangewIC50, './processeddata/dataFoldChange.xlsx', 'Sheet', 'FC significant')
writetable(dataFCpvalue, './processeddata/dataFoldChange.xlsx', 'Sheet', 'p-values')
writetable(table('Source:','load_RPPA_data.m', cd), '../processeddata/dataFoldChange.xlsx', 'Sheet', 'anno')

save mat_dataFoldChange.mat dataFoldChange* dataFC*

