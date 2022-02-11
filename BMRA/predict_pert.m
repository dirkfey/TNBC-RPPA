load mat_bmra_analysis_ConsNW_19-09-09.mat
load mat_BMRAresults_19-09-09.mat 

iCell = 3;
%A = eye(size(myCOnsNW_resi)) + myCOnsNW_resi;
A = BMRAresults(iCell).BS;
A = eye(size(A)) + A;

% dx = Ax + B u => x = A^-1 B u 

u = -1;

B = myNetwork == "AKT_p";
x3 = inv(A)*B

table(x1,x2,x3, 'RowNames', myNetwork)

