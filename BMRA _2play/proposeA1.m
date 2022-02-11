function [A1] = proposeA1(A,A0,NKMAX)
%A1=A;
if(sum(A)==NKMAX)
   A1=remove(A);
else
    if(sum(A)==0)
        A1=add(A,A0);
    else
        toss=rand();
        if(toss<0.5)
          A1=add(A,A0);
        else
           A1=remove(A);
        end
    end
end
    
    function [B] =add(A,A0)
        B=A;  
        i=find((A0-A)==1);
        B(i(ceil(rand()*length(i))))=1;
    end

function [B] =remove(A)
        B=A;  
        i=find(A==1);
        j=ceil(rand()*length(i));
        B(i(j))=0;
    end

end
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% A1=A;
% if sum(A)==min(length(A)-1,n-1)
%     index=find(A==1);
%     j=ceil(rand()*length(index));
%     A1(index(j))=0;
% end
% 
% if sum(A)==0
%     index=setdiff(1:length(A),i);
%     j=ceil(rand()*length(index));
%     A1(index(j))=1;
% end
% 
% if sum(A)>=1 && sum(A)<min(length(A)-1,n-1)
%     toss=rand();
%     if(toss<0.5)
%         index=find(A==1);
%         j=ceil(rand()*length(index));
%         A1(index(j))=0;
%     else
%         index=setdiff(find(A==0),i);
%         j=ceil(rand()*length(index));
%         A1(index(j))=1;
%     end
% end
%     
% end

