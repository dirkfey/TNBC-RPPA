function [ r,p] = sample_r_st(a,bs,p,n,Vs,mus)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
as=a+(n+p)/2;
%as=a+n/2;
sig=(bs/as)*Vs;
%mus
%mvtrnd(sig,as,1);
r1=mvtrnd(sig,as,1);
r=mus'+r1;
p=mvtpdf(r1',sig,as);
end

