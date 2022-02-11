function [T ] = df_lump(T0, lump, lumpedname)

lumpchar = convertStringsToChars(lump);
%keyboard
X = T0{:,lumpchar};
x_lumped = mean(X,2);
T = T0;
%keyboard
T{:,lumpchar{1}} = x_lumped;
T.Properties.VariableNames{T.Properties.VariableNames==lump(1)} = lumpedname;
T(:,lumpchar(2:end)) = [];

%myNetwork(myNetwork==lumpchar(1)) = lumpedname;
%myNetwork(myNetwork==lumpchar(2:end)) = [];

end