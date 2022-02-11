function [lp] = Lp_AGivenY_new(A,A0,Vb_1, Vs, as, bs, delta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[log(det(Vs))+log(det(Vb_1)) log(n/(n+1))]
lp=0.5*log(det(Vs))+0.5*log(det(Vb_1))-as*log(bs)-delta*sum(abs((A0-A)));%+log(nchoosek(na,p))+log(beta(alpha1+p,beta1+na-p));
end

