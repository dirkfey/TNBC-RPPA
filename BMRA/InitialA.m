function [A] = InitialA(A_init,NKMAX)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%A_init
%NKMAX
index=find(A_init==1);
%rand(1,NKMAX)
%length(index)
l=ceil(rand(1,NKMAX)*length(index));
l=unique(l);
A=zeros(size(A_init));
A(index(l))=1;

end

