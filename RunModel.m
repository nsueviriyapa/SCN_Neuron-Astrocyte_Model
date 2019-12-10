clear all; close all; clc; 

rr=1; %# of independent run %%%%%

 
for ii = 1:rr
    

options = odeset('RelTol',1e-3,'AbsTol',1e-6);

ncell_A = 1;
ns_A = 23; 

ncell_N =100;      %# of cells
ns_N = 21;         %# of states

ncell=ncell_N+ncell_A;



%%% Create Heterogeneity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% Neurons with random perturbations %%%

sd = 0.01;   %  percent standard deviation
 
vsP0_n =(0.94*ones(ncell_N,1))+ sqrt((0.94*sd*6)^2)*randn(ncell_N,1);            
vsB_n = (1.0*ones(ncell_N,1))+ sqrt((1.0*sd*1)^2)*randn(ncell_N,1);              
vmB_n = (0.8*ones(ncell_N,1))+ sqrt((0.8*sd*1)^2)*randn(ncell_N,1);             


%%% Astrocyte %%%
vsB_a = 1.0;                                                                     
vmB_a = 0.8;                                                                     


%%% Connectivity Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Perc_VIP =0.2 ;           % Perc_VIP is %VIP producersb=rand(1,ncell);    20% VIPergic cells
bita= 0.05 ;              % probability to add connections for N2N connections

P_B     = 1;              %probability to add connections for VIP N2A
P_B1    = 1;              %probability to add connections for GABA A2N
P_B1_r  = 1;              %probability to add connections for GABA A2N release
P_B1_u  = 1;              %probability to add connections for GABA A2N uptake
P_B2    = 1;              %probability to add connections for Glu A2N


[A ,A1]=adjacency_NN(ncell_N,bita,Perc_VIP);                                                                 %For N2N, A--> VIP and A1-->GABA
[B, B1, B1_r, B1_u, B2] =adjacency_NA_AN(ncell_N,ncell_A,A, P_B, P_B1, P_B1_r, P_B1_u, P_B2);                %For N2A or A2N, B --> VIP


%%% VIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For VIP N2N

sumal=1./sum(A,2)';                      
sumal(isinf(sumal))=0;                   
scale_VIP_N2N = mean(sum(A,2))             % mean number of connections per cells

%For VIP N2A

sumal_VIP_NA=1./sum(sum(B,2));                                     
scale_A2N = mean(sum(B,2))                 % mean number of connections per cells



%%% GABA and Glutamate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For GABA N2N

sumalGABA=1./sum(A1,2)';
sumalGABA(isinf(sumalGABA))=0;
mean(sum(A1,2))
scale_GABA_N2N = mean(sum(A1,2))             % mean number of connections per cells

%For GABA A2N

sumal_GABA_AN_r=1./sum(sum(B1_r,2));  
sumal_GABA_AN_u=1./sum(sum(B1_u,2));  

%For Glu A2N

sumal_Glu_AN=1./sum(sum(B2,2));       


    

%%% Solve ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = [0:0.1:600];
[t,y]=ode23(@ODEs,t,IC1(ncell_N,ncell_A,ncell,ns_N,ns_A),options,Parameters,A,sumal,A1,sumalGABA,ncell_N,ncell_A,ncell,ns_N,ns_A,vsP0_n,vsB_n,vmB_n,B,B1,B1_r, B1_u,B2,sumal_VIP_NA,sumal_Glu_AN,sumal_GABA_AN_r,sumal_GABA_AN_u,vsB_a,vmB_a);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Neuronal and Astrocyic PermRNA outputs 


for j =1:ncell_N
    
     MP_n(:,j)= y(:,(j-1)*ns_N+3);
   
end



for m=1:ncell_A

      MP_a(:,m)= y(:,(ncell_N*ns_N)+(m-1)*ns_A+3);
    
end


end







