function [r,A,LPA,LPR] = sample_corrected_student_t( X,Y,A_init,r_init,a,b,lambda,alpha1,beta1,noit,node )
r=zeros(size(X,2),noit);
A=zeros(size(X,2),noit);
LPA=zeros(noit,1);
LPR=zeros(noit,1);
%index=find(A_init==1);
delta=1;
%%initialize log posterior
NKMAX=min(sum(A_init),size(X,1)-1);
A0=A_init;
A_init=InitialA(A0,NKMAX);
r_init=r_init'.*A_init;
Xe=X(:,logical(A_init));
[c, p, n, V1, Vb_1, Vs_1, Vs, mus, as, bs, na] = Prepare(Xe,Y,A_init,a,b,lambda);
%Lp_A_Y=Lp_AGivenY(p, n,  Vb_1, Vs, as, bs,alpha1,beta1, na);
Lp_A_Y=Lp_AGivenY_new(A_init, A0,  Vb_1, Vs, as, bs,delta);
for i=1:noit
    temp_A=proposeA1(A_init,A0,NKMAX);
    Xe_temp=X(:,logical(temp_A));
    [c_t, p_t, n_t, V1_t, Vb_1_t, Vs_1_t, Vs_t, mus_t, as_t, bs_t, na_t] = Prepare(Xe_temp,Y,temp_A,a,b,lambda);
    Lp_A_Y_t=Lp_AGivenY_new(temp_A,A0,Vb_1_t, Vs_t, as_t, bs_t,delta);
    if(Lp_A_Y_t>Lp_A_Y)
        A_init=temp_A;
        Lp_A_Y=Lp_A_Y_t;
        c=c_t;
        p=p_t;
        n=n_t;
        V1=V1_t;
        Vb_1=Vb_1_t;
        Vs_1=Vs_1_t;
        Vs=Vs_t;
        mus=mus_t;
        as=as_t;
        bs=bs_t;
        na=na_t;
    else
        if sample_mini([Lp_A_Y_t Lp_A_Y],1)==1
            A_init=temp_A;
            Lp_A_Y=Lp_A_Y_t;
            c=c_t;
            p=p_t;
            n=n_t;
            V1=V1_t;
            Vb_1=Vb_1_t;
            Vs_1=Vs_1_t;
            Vs=Vs_t;
            mus=mus_t;
            as=as_t;
            bs=bs_t;
            na=na_t;
        end
    end

    A(:,i)=A_init;
    r1=zeros(size(r_init));
    %sum(A_init)
    %sample_r_st(a,bs,p,n,Vs,mus)
    [rt,rp]=sample_r_st(a,bs,p,n,Vs,mus);
    r1(logical(A_init))=rt;
    LPA(i)=Lp_A_Y;
    LPR(i)=rp;
    r(:,i)=r1;
end

end

