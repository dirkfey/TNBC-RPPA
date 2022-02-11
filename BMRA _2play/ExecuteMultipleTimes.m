function [PS,BS,CS,Rs,AS] = ExecuteMultipleTimes( data, pert,G,noit,burnin, times)
%% Inputs:  
%%(a)data= (N x Np) Global response matrix; N= Number of nodes; Np = Number of perturbation experiments; 
%%(b)pert = N x Np binary matrix where pert(i,j)=1 indicates that data(:,j) is the global response of the network when node i node is perturbed.
%%(c)noit = number of Gibbs scans; 
%%(d)times =number of times the program is intended to be execute
%% Outputs:
%%(a) PS= posterior edge probability matrix .


a=1;
b=1;
lambda=0.2;
alpha1=1;
beta1=1;

sd=size(data);
AS=zeros(sd(1),sd(1),(noit-burnin)*times);
RS=zeros(sd(1),sd(1),(noit-burnin)*times);


PS=zeros(sd(1),sd(1));
BS=zeros(sd(1),sd(1));
CS=zeros(sd(1),sd(1));


for i=1:times
    [As,Rs,LPAs,LPRs]=Net_inf_corrected(data,pert,G,a,b,lambda,alpha1,beta1,noit);
    AS(:,:,(i-1)*(noit-burnin)+1:i*(noit-burnin))=As(:,:,burnin+1:noit);
    RS(:,:,(i-1)*(noit-burnin)+1:i*(noit-burnin))=Rs(:,:,burnin+1:noit);

end
PS=mean(AS,3);
%for j=1:sd(1)
for i=1:sd(1)
    for j=1:sd(1)
       r=reshape(RS(i,j,:),1,size(RS,3));
       %DF:
       r1 = r;
       %r1=outliers(r,size(r,3)/10);
       noutliers = ceil(length(r1)/20);
       [~, sidx] = sort(r1); r1(sidx([1:noutliers (end-noutliers):end])) = [];
       %r1(r1==0)=[];
       %r1(isoutlier(r1))=[];
       %r1=r1(~isnan(r1));
       BS(i,j)=mean(r1);%meanr(RS(j,:,:),LPAs(j,:),LPRs(j,:));
       CS(i,j)=std(r1);%stdr(RS(j,:,:),LPAs(j,:),LPRs(j,:));
       %RS(j,(i-1)*(noit-burnin)+1:i*(noit-burnin),:)=Rs(j,burnin+1:noit,:);
    end
end
%PS=PS/times;
%BS=BS/times;
%CS=CS/times;
%Rs=reshape(RS,size(RS,1)*size(RS,2),size(RS,3))';
Rs = RS;
end

