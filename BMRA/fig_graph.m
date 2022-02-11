
%These are the interaction in the svg file:
% Format: inkscape:label  from-node  to-node

mysvgInteractions = [
    "RAF-MEK"   "cRAF_S338_"       "MEK1_2_S217_" 
    "MEK-ERK"   "MEK1_2_S217_"     "MAPK_T202_"   
    "ERK-RAF"   "MAPK_T202_"       "cRAF_S338_"         
    "ERK-PTEN"  "MAPK_T202_"       "PTEN"         
    "PDK-AKT"   "PDK1_S241_"       "AKT_p" 
    "AKT-PTEN"  "AKT_p"            "PTEN"         
    "MTOR-RAF"  "mTOR_p"           "cRAF_S338_"   
    "MTOR-S6"   "mTOR_p"           "S6RIB_p"      
    "MTOR-NFkB" "mTOR_p"           "NF_kB_p65"    
    "MTOR-CHK1" "mTOR_p"           "CHK1_S345_"   
    "MTOR-p38"  "mTOR_p"           "p38MAPK_T180_"
    "MTOR-FAK"  "mTOR_p"           "FAK_Y925_"    
    "GSK3B-PTEN"    "GSK3B_S9_"        "PTEN"         
    "GSK3B-NFkB"    "GSK3B_S9_"        "NF_kB_p65"    
    "AMPK-FAK"  "AMPK_T172_"       "FAK_Y925_"    
    "NFkB-PTEN" "NF_kB_p65"        "PTEN"         
    "NFkB-MTOR" "NF_kB_p65"        "mTOR_p"       
    "p38-PTEN"  "p38MAPK_T180_"    "PTEN"         
    "STAT3-MEK" "stat3_y750_"      "MEK1_2_S217_" 
    "STAT3-ERK" "stat3_y750_"      "MAPK_T202_"   
    "STAT3-PTEN"    "stat3_y750_"      "PTEN"         
    "STAT3-FAK" "stat3_y750_"      "FAK_Y925_" 
    "PI3K-PDK" "PI3KalPha"        "PDK1_S241_"   
    "GSK3B-S6"  "GSK3B_S9_"        "S6RIB_p"      
    "S6K-S6"    "P70_S6K_T389_"    "S6RIB_p"      
    "S6K-p38"   "P70_S6K_T389_"    "p38MAPK_T180_"
    "NFkB-AMPK" "NF_kB_p65"        "AMPK_T172_"   
    "NFkB-p53"  "NF_kB_p65"        "p53"          
    "NFkB-FAK"  "NF_kB_p65"        "FAK_Y925_"    
    "STAT3-p38" "stat3_y750_"      "p38MAPK_T180_"
    "STAT3-p53" "stat3_y750_"      "p53"          
    "FAK-ERK"   "FAK_Y925_"        "MAPK_T202_"   
    "FAK-NFkB"  "FAK_Y925_"        "NF_kB_p65"    
    "FAK-p53"   "FAK_Y925_"        "p53" 
    "MTOR-S6K"  "mTOR_p"           "P70_S6K_T389_"
    "S6K-MTOR"  "P70_S6K_T389_"    "mTOR_p"          
    "AMPK-MTOR" "AMPK_T172_"       "mTOR_p"
    "CHK1-p53"  "CHK1_S345_"       "p53"
    "p38-p53"   "p38MAPK_T180_"    "p53"
    "MTOR-AKT"  "mTOR_p"           "AKT_p" 
    "AKT-MTOR"   "AKT_p"            "mTOR_p"            
    "AKT-GSK3B" "AKT_p"            "GSK3B_S9_"        
    ];

myNW = myInteractionMatrix;
xDoc = xmlread('./graph_NW_2_base.svg');
xDoc = svg_colorInteractions(xDoc, mysvgInteractions, myNW, myNetwork, [0 0 0]);
xmlwrite('../figs/graph_NW_2_canonical_tmp.svg',xDoc)


myNW = myConsNW_sens;
myNW(abs(myNW)<0.4)=0;
xDoc = xmlread('./graph_NW_2_base.svg');
xDoc = svg_colorInteractions(xDoc, mysvgInteractions, myNW, myNetwork);
xmlwrite('../figs/graph_NW_2_sensitive_tmp.svg',xDoc)

myNW = myCOnsNW_resi;
myNW(abs(myNW)<0.4)=0;
xDoc = xmlread('./graph_NW_2_base.svg');
xDoc = svg_colorInteractions(xDoc, mysvgInteractions, myNW, myNetwork);
xmlwrite('../figs/graph_NW_2_resistant_tmp.svg',xDoc)


function xDoc = svg_colorInteractions(xDoc, mysvgInteractions, myNW, myNetwork, varargin )

idx_interactions = find(myNW~=0);
[i_to, i_from] = ind2sub(size(myNW), idx_interactions);
to = myNetwork(i_to);
from = myNetwork(i_from);

%loop through all interactions
svgcolor = [0 0 0];
myCM = lines(4);
for i=1:length(i_to)
    maskto = mysvgInteractions(:,3)==to(i);
    maskfrom = mysvgInteractions(:,2)==from(i);
    if myNW(i_to(i),i_from(i)) > 0 
        svgcolor = myCM(2,:);
    else
        svgcolor = myCM(1,:);
    end
    if ~isempty(varargin) % overrides above sign based colors
        svgcolor = varargin{1};
    end
    idx = find(maskto&maskfrom);
    if ~isempty(idx)
        svglabel = mysvgInteractions(idx,1);
        fprintf('%d %s\n',i,join(["Interaction to" to(i) "from" from(i) "found"]))
        %set color:
        xDoc = svg_setPathColorByLabel(xDoc, svglabel,svgcolor);
    else
        fprintf(join(["   Interaction to" to(i) "from" from(i) "not found\n"]))
    end
end

end
