function [A A1]=adjacency_NN(ncell_N,bita,Perc_VIP)

kk=1:ncell_N;  
X=kk(ones(1,ncell_N),:); 
s=sqrt(ncell_N); 
kk=kk';

Var1=fix((kk-1)/s);  
Var1=Var1(:,ones(1,ncell_N));

Var2=rem((kk-1),s);  
Var2=Var2(:,ones(1,ncell_N));

A = (sqrt((Var1-fix((X-1)/s)).^2+(Var2-rem((X-1),s)).^2)).\1;
A(A~=1)=0;

clear Var1 Var2 

%CREATE PERIODIC BOUNDARY CONDITIONS to eliminate end effects at the edges of the lattice

%FOLD THE TOPS

for i=1:s   
A(i,s*(s-1)+i)=1;
A(s*(s-1)+i, i)=1;
end

%FOLD SIDES  %Create Vector with Specified Increment

w=1:s:ncell_N; 
for o=1:length(w)
A(w(o),(s-1)+w(o))=1;
A((s-1)+w(o),w(o))=1;
end

 
%EXCLUDE SELF CONNECTIONS

for i=1:ncell_N
    A(i,i)=NaN; 
end


%ADD CONNECTIONS RANDOMLY

for i=1:ncell_N
    r=X(i,A(i,:)==0);     
    b=rand(length(r),1);  
    b(b<bita)=0;
    r=r(b==0);
    if ~isempty(r)
        for k=1:length(r)
            A(i,r(k))=1;  
            A(r(k),i)=1; 
        end
    end
end

A(isnan(A))=0; % REPLACE NaN WITH 0

A1=A;
clear b r k i w X 


%VIP PERCENTAGE (VIPergic Cells/Erase Connections from NON-VIPergic)

b=rand(1,ncell_N);    
b(b <Perc_VIP)=0;      
b(b~=0)=1;
wa=kk(b==1);

    for o = 1:length(wa)
        A(:,wa(o))=0;  
    end
    
clear Perc_VIP b wa s kk
    
   
