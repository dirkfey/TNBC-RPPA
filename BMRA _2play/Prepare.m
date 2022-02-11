function [c, p, n, V1, Vb_1, Vs_1, Vs, mus, as, bs, na] = Prepare(X,Y,A,a,b,lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c=length(Y);
p=size(X,2);
n=length(Y);
V1=X'*X;
Vb_1=(V1+lambda*eye(size(X,2)))/(c);
Vs_1=V1+Vb_1;
Vs=Vs_1\eye(size(Vs_1));
%size(Vs)
%size(X)
%size(Y)
%lambda
mus=Vs*X'*Y;
as=a+n/2;
bs=b+0.5*(Y'*Y-mus'*Vs_1*mus);

na=length(A)-1;

end

