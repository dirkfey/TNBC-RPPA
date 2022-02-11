%bmra_prepare_data ... Skript to set up BMRA
%   list of drugs, list of analytes, network, network interactions, etc

%list of Analytes:
myAnalytes = [
    "cRAF_S338_"
    'MEK1_2_S217_'
    'MAPK_T202_'
    'PI3KalPha'
    'PTEN'
    'PDK1_S241_'
    'AKT_S473_' 
    'AKT_T308_'
    'mtor_s2481_'
    'GSK3B_S9_'
    'AMPK_T172_'
    'mTOR_S2448_'
    'P70_S6K_T389_'
    'S6RIB_S235_'
    'S6RIB_S240_'
    'NF_kB_p65'
    'CHK1_S345_'
    'p38MAPK_T180_'
    'p53'
    'stat3_y750_'
    'FAK_Y925_'
    ];

%List of network nodes (some analytes will be lumped into one network node!)
myNetwork = [
    "cRAF_S338_"
    'MEK1_2_S217_'
    'MAPK_T202_'
    'PI3KalPha'
    'PTEN'
    'PDK1_S241_'
    'AKT_p' 
    'mTOR_p'
    'GSK3B_S9_'
    'AMPK_T172_'
    'P70_S6K_T389_'
    'S6RIB_p'
    'NF_kB_p65'
    'CHK1_S345_'
    'p38MAPK_T180_'
    'p53'
    'stat3_y750_'
    'FAK_Y925_'
    ];
nNodes = length(myNetwork);

%List of canonical interactions and corresponding interaction matrix:
myInteractionMatrix = zeros(nNodes);
myInteractionsList = [
    1 3     %RAF<-ERK
    2 1     %MEK<-RAF
    3 2     %ERK<-MEK
    5 15    %PTEN<-NFkB
    7 5     %AKT<-PTEN
    6 4     %PDK<-PI3K
    7 6     %AKT-T308<-PDK
    7 8     %AKT-S473<-mTORC2
    8 11    %mTORC2<-S6K
    9 7    %GSK3B<-AKT-S473
    8 10   %mTORC1<-AMPK
    %8 9    %mTORC1<-GSK3B
    8 7    %mTORC1<-AKT-S473
    11 8   %S6K<-mTORC1
    12 11   %S6RIB<-S6K
    15 8   %NFkB<-mTORC1
    16 14   %p53<-CHK1   
    16 15   %p53<-p38
    ];
linearIdx = sub2ind([nNodes nNodes],myInteractionsList(:,1),myInteractionsList(:,2));
myInteractionMatrix(linearIdx)=1;

%Mapping of drugs to target (network node)
myDrugs = [
    ... %'AICAR'       "AMPK_T172_"
    'CHIR-98014'  "GSK3B_S9_"
    'Dorsomorphin' "AMPK_T172_"
    'Everolimus'  "mTOR_p"
    'Ipatasertib' "AKT_p"
    'PF-00562271' "FAK_Y925_"
    'PF-4708671'  "P70_S6K_T389_"
    'QNZ'         "NF_kB_p65"
    'Stattic'     "stat3_y750_"
    'TAK-632'     "cRAF_S338_"
    'U0126'       "MAPK_T202_" 
    ]; 

[~, myPertubedNodes] = ismember(myDrugs(:,2), myNetwork);

nNodes = length(myNetwork);
nPert = size(myDrugs,1);

%perturbation matrix needed for BMRA
myPert = zeros(nNodes, nPert);
for ipert = 1:nPert
    inode = find(myNetwork==myDrugs(ipert,2));
    myPert(inode, ipert) = 1;
end

clear inode ipert linearIdx 
