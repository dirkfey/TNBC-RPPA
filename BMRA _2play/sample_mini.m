function [ i] = sample_mini( l,T )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% p0=exp(l(1));
% p1=exp(l(2));
% sum=p0+p1;
% P0=[p0/sum p1/sum];
l=l/T;
f=max(l);
logsum=f+log(sum(exp(l-f)));
P=cumsum(exp(l-logsum));
r=rand();
i=min(find(P>=r));
end

