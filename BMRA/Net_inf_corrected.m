function [As,Rs,LPAs,LPRs ] = Net_inf_corrected(data,pert,G,a,b,lambda,alpha1,beta1,noit)
%data=Number of measured nodes X number of perturbations 
%pert=perturbation matrix
Nodes=size(data,1);
As=zeros(Nodes,Nodes,noit);
Rs=zeros(Nodes,Nodes,noit);
LPAs=zeros(Nodes,noit);
LPRs=zeros(Nodes,noit);
%mrb=zeros(Nodes);
%srb=zeros(Nodes);
%G=createInitialGraph1(size(data,1),size(data,2)-1);

rr=randn(size(data,1)).*G;
%try
parfor i=1:Nodes
%for i=1:Nodes
    Ai=G(i,:);
    if(sum(Ai)>0)
        
    ri=rr(i,:)';
    XY=data(:,logical(1-pert(i,:)));
    Y=XY(i,:)';
    X=XY';
    %initialize Ai and ri
    %make sure the number of nonzero connection coefficients does not
    %exceed the number of perturbation experiments
    if(sum(Ai)>length(Y))
        %Tapesh:
        %indx=find(ri==1);
        %indxe=ceil(rand(1,sum(ri))-length(Y))*length(indx));
        %Dirk:
        indx=find(~(ri==0));
        indxe=ceil(rand(1,sum(~(ri==0))-length(Y))*length(indx));
        %
        ri(indx(indxe))=0;
        Ai(indx(indxe))=0;
    end
   % Ai
    %start sampling
    [r, A, LPA, LPR]=sample_corrected_student_t( X,Y,Ai,ri,a,b,lambda,alpha1,beta1,noit,i);
    %store data
    %mr=mean(r(burnin+1:noit,:));
    %sr=std(r(burnin+1:noit,:));
    %pi=mean(A(burnin+1:noit,:));
    %P(i,:)=pi;
    %mrb(i,:)=mr;
    %srb(i,:)=sr;
    As(i,:,:)=A;
    Rs(i,:,:)=r;
    LPAs(i,:)=LPA;
    LPRs(i,:)=LPR;
    end
end
%catch
%    keyboard
%end

end

