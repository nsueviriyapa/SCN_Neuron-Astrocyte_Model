function IC2 = IC1(ncell_N,ncell_A,ncell,ns_N,ns_A)

z=(ns_N*ncell_N)+(ns_A*ncell_A);
i=1:ncell;
IC2=zeros(1,z); 

%% Neurons

for i = 1:ncell_N
    
        IC1(1) =0.1;
        IC1(2)= 0.1;
        IC1(3) = 2.79566;
        IC1(4) = 1.964254;
        IC1(5) = 7.946372;
        IC1(6) = 0.4;
        IC1(7) = 12;
        IC1(8) = 0.126363;
        IC1(9) = 8.986265;
        IC1(10) = 1.25558;
        IC1(11) = 0.165471;
        IC1(12) = 0.197619;
        IC1(13) = 0.090808;
        IC1(14) = 2.408975;
        IC1(15) = 0.479535;
        IC1(16) = 1.941518;
        IC1(17) = 0.32736;
        IC1(18) = 0.049602;
        IC1(19)= 0.12;
        IC1(20)=0;
        IC1(21)=0.01;
         
    
    for j = 1:ns_N 
        IC2((i-1)*ns_N+j) = IC1(j);
    end
    
end


%% Astrocyte


for i=ncell_N+1:ncell
    
        IC1(1) =0.1;
        IC1(2)= 0.1;
        IC1(3) = 2.79566;
        IC1(4) = 1.964254;
        IC1(5) = 7.946372;
        IC1(6) = 0.4;
        IC1(7) = 12;
        IC1(8) = 0.126363;
        IC1(9) = 8.986265;
        IC1(10) = 1.25558;
        IC1(11) = 0.165471;
        IC1(12) = 0.197619;
        IC1(13) = 0.090808;
        IC1(14) = 2.408975;
        IC1(15) = 0.479535;
        IC1(16) = 1.941518;
        IC1(17) = 0.32736;
        IC1(18) = 0.049602;
        IC1(19)= 0.12;
        IC1(20)=0.00;  
        IC1(21)=0.00;
        IC1(22)=0.00;
        IC1(23)=0.00;
    
   
        for j = 1:ns_A 
            IC2((ncell_N*ns_N)+((i-ncell_N-1)*ns_A)+j) = IC1(j);
        end
    
end




