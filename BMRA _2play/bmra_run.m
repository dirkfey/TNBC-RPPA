% bmra_run ... Skript to run BMRA analysis on all cell lines
%
%   Input: 
%   Outputs: 
%   Dependencies: df_lump.m, ExecuteMultipleTimes.m, Net_inf_corrected.m,
%       sample_corrected_student_t.m, InitialA.m, Prepare.m,
%       Lp_AGivenY_new.m, proposeA1.m, sample_mini.m, sample_r_st.m, 
%       df_pcolor.m

hyperp.noit = 1000%10000;
hyperp.burnin = 500;
hyperp.times = 10^2;%10*4;

load ../mat_dataFoldChange.mat
bmra_prepare_data

myCellLines = ["MDA-MB-468" "CAL-85-1" "MDA-MB-231" 'MCF10A' 'HDQP-1' 'MDA-MB-157' 'MFM223'];    

%myAnalytes = myNetwork';
   
clear BMRAresults
figure(99), clf
for i=1:length(myCellLines)
    % select cell line:
    myCellLine = myCellLines(i);
    %mySelection0 = dataFoldChangewIC50(dataFoldChangewIC50.CellLine==myCellLine,:);
    mySelection0 = dataFoldChange_allwIC50(dataFoldChangewIC50.CellLine==myCellLine,:);
    % log2 transform fold changes
    mySelection0{:,3:end-1} = log2(mySelection0{:,3:end-1});

%Impute INFs
tmp = (mySelection0{:,3:(end-1)});
tmp(isinf(tmp)) = max(max(tmp(~isinf(tmp))));
tmp(isnan(tmp)) = 1;
mySelection0{:,3:(end-1)} = tmp

[mask, idxAnalytes] = ismember(['CellLine'; 'Treatment'; convertStringsToChars(myAnalytes(:)); 'IC50'], mySelection0.Properties.VariableNames);
mySelection0 = mySelection0(:,idxAnalytes(idxAnalytes>0));
mask = ismember(mySelection0.Treatment,convertStringsToChars(myDrugs(:,1)));
mySelection0(~mask,:) = [];
[mySelection0] = df_lump(mySelection0, ["AKT_S473_" "AKT_T308_"], 'AKT_p');
[mySelection0] = df_lump(mySelection0, ["mtor_s2481_" "mTOR_S2448_"], 'mTOR_p');
[mySelection0] = df_lump(mySelection0, ["S6RIB_S235_" "S6RIB_S240_"], 'S6RIB_p');

%Convert to DataMatrix 
import bioma.data.*
mySelection00 = sortrows(mySelection0, 'IC50');
myDM = DataMatrix(mySelection00{:,3:end},'RowNames', mySelection00.Treatment,'ColNames', mySelection00.Properties.VariableNames(3:end) );

clear PS BS CS Rs 
%filter high values out (mostly PTEN):
%myDMM= double(myDM);
%myDMM(myDMM>3) = 0;
%myDMM(myDMM<-3) = 0;
%myDM(:,:) = myDMM;

%Set up global response matrix 
R = [];
for ii=1:nNodes;
    R = [R mySelection00{:,convertStringsToChars(myNetwork(ii))}];
end
%R(R>3)  = 0;
%R(R<-3) = 0;
R = R';

%set up interaction matrix 
G0 = myInteractionMatrix;
%G0(G0==0)=1; all interactions 
for ii = 1:nNodes
    G0(ii,myPertubedNodes) = 1; %all interactions from perturbed nodes 
end
G0(:,1) = myInteractionMatrix(:,1); %Raf phos only MEK 
G0(:,2) = myInteractionMatrix(:,2); %MEK only ERK
G0(diag(ones(nNodes,1))==1)=0;

readme = '[PS,BS,CS,Rs] = ExecuteMultipleTimes(R, myPert, G0, hyperp.noit , hyperp.burnin, hyperp.times);';
[PS,BS,CS,Rs] = ExecuteMultipleTimes(R, myPert, G0, hyperp.noit , hyperp.burnin, hyperp.times);

clear Rs

BMRAresults(i).CellLine=myCellLine;
BMRAresults(i).PS=PS;
BMRAresults(i).BS=BS;
BMRAresults(i).CS=CS;
BMRAresults(i).R = R;
BMRAresults(i).data = mySelection00;
BMRAresults(i).DM = myDM;

%
figure(99), 
df_pcolor(BS.*(PS>0.0))
hold on
[ix, iy] = find(myInteractionMatrix'==1);
scatter(ix+0.5, iy+0.5, 'xk')
colormap(redbluecmap)
set(gca, 'xTick', [1:nNodes]+0, 'XTickLabel', myNetwork, 'XTickLabelRotation', 60)
set(gca, 'yTick', [1:nNodes]+0.5, 'yTickLabel', myNetwork)
set(gca, 'CLim', [-1 1])
set(gca, 'YDir', 'reverse')

end

save mat_BMRAresults.mat BMRAresults readme myPert G0 hyperp myCellLines
    

%% Testing on interactions:
myIntMat_sens = [
    reshape(BMRAresults(1).BS',1,[])
    reshape(BMRAresults(2).BS', 1, [])
    reshape(BMRAresults(3).BS',1,[])
    ];
myIntMat_resi = [
    reshape(BMRAresults(5).BS',1,[])
    reshape(BMRAresults(6).BS',1,[])
    reshape(BMRAresults(7).BS',1,[])
    ];
tmp = repmat(myNetwork, 1, nNodes);
myIntMat_labels1 = reshape(tmp,1,[])
myIntMat_labels2 = reshape(tmp',1,[])
myIntMat_labels  = join([myIntMat_labels1', myIntMat_labels2'], ' -> ')

%filter all zeros
mask = all(myIntMat_sens==0);
myIntMat_sens(:,mask)=[];
myIntMat_resi(:,mask)=[];
myIntMat_labels(mask)=[];
clear P D
for i=1:size(myIntMat_sens,2)
    [~, P(i)] = ttest2(myIntMat_sens(:,i),myIntMat_resi(:,i));
    %P(i) = ranksum(myIntMat_sens(:,i),myIntMat_resi(:,i));
    D(i) = mean(myIntMat_resi(:,i))-mean(myIntMat_sens(:,i));
end

myTest = table(myIntMat_labels(:), P(:), D(:), 'VariableNames', {'Interaction' 'pvalue' 'Delta'});
myTest = [myTest ... 
    array2table([myIntMat_sens(:,:)' myIntMat_resi(:,:)'], 'VariableNames', convertStringsToChars(replace(myCellLines([1 2 3 5 6 7]), "-", "_")) ) ];
myTest = sortrows(myTest, 'pvalue', 'ascend')

writetable(myTest,'NW_changes_sensVSresi.xlsx', 'Sheet', 'NW interactions')
writetable(table("Source:", "/code/BMRA_2play/bmra_run.m"),'NW_changes_sensVSresi.xlsx', 'Sheet', 'anno')
fileID = fopen('NW_changes_sensVSresi_readme.txt','w');
fprintf(fileID, 'Source: /code/BMRA_2play/bmra_run.m\n');
fprintf(fileID, 'Date: %s\n', datestr(now));
fclose(fileID);

figure(6), clf
scatter(D,-log10(P))

hold on
plot(8*[-1 1], -log10([0.05 0.05]), 'k--')

