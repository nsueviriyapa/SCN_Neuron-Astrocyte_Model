function dydt= ODEs(t,y,p,A,sumal,A1,sumalGABA,ncell_N,ncell_A,ncell,ns_N,ns_A,vsP0_n,vsB_n,vmB_n,B,B1,B1_r ,B1_u,B2,sumal_VIP_NA,sumal_Glu_AN,sumal_GABA_AN_r,sumal_GABA_AN_u,vsB_a,vmB_a)

  z=(ns_N*ncell_N)+(ns_A*ncell_A);

  dydt=zeros(z,1); 

  j=1:ncell_N;
  m=1:ncell_A;


%%% Neuronal cells

Ca_n=y((j-1)*ns_N+1);
Ca_store_n =  y((j-1)*ns_N+2);
MP_n = y((j-1)*ns_N+3);
MC_n = y((j-1)*ns_N+4);
MB_n = y((j-1)*ns_N+5);
PC_n = y((j-1)*ns_N+6);
CC_n = y((j-1)*ns_N+7);
PCP_n = y((j-1)*ns_N+8);
CCP_n = y((j-1)*ns_N+9);
PCC_n = y((j-1)*ns_N+10);
PCN_n = y((j-1)*ns_N+11);
PCCP_n = y((j-1)*ns_N+12);
PCNP_n = y((j-1)*ns_N+13);
BC_n = y((j-1)*ns_N+14);
BCP_n = y((j-1)*ns_N+15);
BN_n = y((j-1)*ns_N+16);
BNP_n = y((j-1)*ns_N+17);
IN_n = y((j-1)*ns_N+18);
CB_n = y((j-1)*ns_N+19);
vVIP_n = y((j-1)*ns_N+20);
gGABA_n=y((j-1)*ns_N+21);


%%% Astrocytic cells

Ca_a=y((ncell_N*ns_N)+(m-1)*ns_A+1);
Ca_store_a = y((ncell_N*ns_N)+(m-1)*ns_A+2);
MP_a = y((ncell_N*ns_N)+(m-1)*ns_A+3);
MC_a = y((ncell_N*ns_N)+(m-1)*ns_A+4);
MB_a = y((ncell_N*ns_N)+(m-1)*ns_A+5);
PC_a = y((ncell_N*ns_N)+(m-1)*ns_A+6);
CC_a = y((ncell_N*ns_N)+(m-1)*ns_A+7);
PCP_a = y((ncell_N*ns_N)+(m-1)*ns_A+8);
CCP_a = y((ncell_N*ns_N)+(m-1)*ns_A+9);
PCC_a = y((ncell_N*ns_N)+(m-1)*ns_A+10);
PCN_a = y((ncell_N*ns_N)+(m-1)*ns_A+11);
PCCP_a = y((ncell_N*ns_N)+(m-1)*ns_A+12);
PCNP_a = y((ncell_N*ns_N)+(m-1)*ns_A+13);
BC_a = y((ncell_N*ns_N)+(m-1)*ns_A+14);
BCP_a = y((ncell_N*ns_N)+(m-1)*ns_A+15);
BN_a = y((ncell_N*ns_N)+(m-1)*ns_A+16);
BNP_a = y((ncell_N*ns_N)+(m-1)*ns_A+17);
IN_a = y((ncell_N*ns_N)+(m-1)*ns_A+18);
CB_a = y((ncell_N*ns_N)+(m-1)*ns_A+19);
gGABA_a=y((ncell_N*ns_N)+(m-1)*ns_A+20);
gGABA_store_a=y((ncell_N*ns_N)+(m-1)*ns_A+21);
gGlu_a=y((ncell_N*ns_N)+(m-1)*ns_A+22);
gGlu_store_a=y((ncell_N*ns_N)+(m-1)*ns_A+23);


%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%


%%% Neuron 

k1_n = p(1);
k2_n = p(2);
k3_n = p(3);
k4_n = p(4);
k5_n = p(5);
k6_n = p(6);
k7_n = p(7);
k8_n = p(8);
KAP_n = p(9);
KAC_n = p(10);
KIB_n = p(11);
kdmb_n = p(12);
kdmc_n = p(13);
kdmp_n = p(14);
kdnc_n = p(15);
kdn_n = p(16);
Kd_n = p(17);
Kdp_n = p(18);
Kp_n = p(19);
KmB_n = p(20);
KmC_n = p(21);
KmP_n = p(22);
ksB_n = p(23);
ksC_n = p(24);
ksP_n = p(25);
n_n = p(26);
m_n = p(27);
Vphos_n = p(28);
V1B_n = p(29);
V1C_n = p(30); 
V1P_n = Vphos_n; 
V1PC_n = Vphos_n; 
V2B_n = p(33);
V2C_n = p(34);
V2P_n = p(35);
V2PC_n = p(36);
V3B_n = p(37);
V3PC_n = Vphos_n;
V4B_n = p(39);
V4PC_n = p(40);
vdBC_n = p(41);
vdBN_n = p(42);
vdCC_n = p(43);
vdIN_n = p(44);
vdPC_n = p(45);
vdPCC_n = p(46);
vdPCN_n = p(47);
%vmB_n= p(48);
vmC_n = p(49);
vmP_n = p(50);
%vsB_n = p(51);
vsC_n = p(52);
%vsP_n=p(53);
%RT_n=p(54);
KD_n = p(55);
%k_n=p(56);
%v0_n=p(57);
%v1_n=p(58);
vP_n = p(59);
VMK_n=p(60);
% Ka_n=p(61);
K_1_n = p(62);
K_2_n = p(63);
WT_n = p(64);
CT_n = p(65);
KC_n = p(66);
%vsP0_n = p(67);
%v2_n=p(68);
kf_n=p(69);
IP3_n=p(70);
VM3_n= p(71);
M3_n= p(72);
KR_n= p(73);
KA_n = p(74);
pA_n = p(75);
VM2_n=p(76);
K2_n=p(77);
M2_n= p(78);
kMK_n=p(79);
V_b_n=p(80);
k_b_n=p(81);

%%% Astrocyte

k1_a = p(82);
k2_a = p(83);
k3_a = p(84);
k4_a = p(85);
k5_a = p(86);
k6_a = p(87);
k7_a = p(88);
k8_a = p(89);
KAP_a = p(90);
KAC_a = p(91);
KIB_a = p(92);
kdmb_a = p(93);
kdmc_a = p(94);
kdmp_a = p(95);
kdnc_a = p(96);
kdn_a = p(97);
Kd_a = p(98);
Kdp_a = p(99);
Kp_a = p(100);
KmB_a = p(101);
KmC_a = p(102);
KmP_a = p(103);
ksB_a = p(104);
ksC_a = p(105);
ksP_a = p(106);
n_a = p(107);
m_a = p(108);
Vphos_a = p(109);
V1B_a = p(110);
V1C_a = p(111);
V1P_a = Vphos_a; 
V1PC_a = Vphos_a; 
V2B_a = p(114);   
V2C_a = p(115);
V2P_a = p(116);
V2PC_a = p(117);
V3B_a = p(118);
V3PC_a = Vphos_a;  
V4B_a = p(120);
V4PC_a = p(121);
vdBC_a = p(122);
vdBN_a = p(123);
vdCC_a = p(124);
vdIN_a = p(125);
vdPC_a = p(126);
vdPCC_a = p(127);
vdPCN_a = p(128);
%vmB_a= p(129);
vmC_a = p(130);
vmP_a = p(131);
%vsB_a = p(132);
vsC_a = p(133);
%vsP_a=p(134)
%RT_a=p(135)
KD_a = p(136);
%k_a=p(137)
%v0_a =p(138)
%v1_a=p(139)
vP_a = p(140);
VMK_a=p(141);
%Ka_a=p(142)
K_1_a = p(143);
K_2_a = p(144);
WT_a = p(145);
CT_a = p(146);
KC_a = p(147);
vsP0_a = p(148);
%v2_a = p(149);
kf_a=p(150);
IP3_a=p(151);
VM3_a= p(152);
M3_a= p(153);
KR_a= p(154);
KA_a = p(155);
pA_a = p(156);
VM2_a=p(157);
K2_a=p(158);
M2_a= p(159);
kMK_a=p(160);
V_b_a=p(161);
k_b_a=p(162);

%addtional parameters

v_vo_n = p(163);
K_vo_n = p(164);
n_vo_n = p(165);

v_kk_n = p(166);
K_kk_n = p(167);
n_kk_n = p(168);


v_vo_a = p(169);
K_vo_a = p(170);
n_vo_a = p(171);

v_kk_a = p(172);
K_kk_a = p(173);
n_kk_a = p(174);

nGluR = p(175);
Km_GluR = p(176);
vm_GluR = p(177);

vo_GABA_a = 1*p(178);         

vm_GABA_a = 1*p(179);
n_GABA_a =  p(180);   
Km_GABA_a = p(181);         

k_dGABAa = p(182);     
n_dGABAa = p(183);

k_dGABAa_store = p(184);          
n_dGABAa_store = p(185);

vm_iGAT =  1*p(186); 
n_GAT_a =  p(187);       
Km_GAT_a = p(188);        
n_iGAT =   p(189);        
Km_iGAT =  p(190);

vo_Glu_a = p(191);
vm_Glu_a = p(192);
n_Glu_a =  p(193);
Km_Glu_a = p(194);
 
 
k_dGlua =  p(195);
n_dGlua =  p(196);
 
k_dGlua_store = p(197);
n_dGlua_store = p(198);
 

vm_iEAAT = p(199); 
n_iEAAT = p(200);       
Km_iEAAT = p(201);              
 
n_EAAT_a = p(202);       
Km_EAAT_a = p(203);


%% Neurotransmitters %%%%%%%%%%%%%%%%%%%%%%%%%%


%Neurons


%no neurotransmitter release and coupling for the first 150h 

    if t<150
      switz=0;
      beta1=zeros(ncell_N,1);
      S_GABA=0.1*ones(ncell_N,1);   
      S_GABAp=0.1*ones(ncell_N,1);  
      beta2=zeros(ncell_N,1);
     
      
    else 
        switz=1;

    %%% VIP calcs %%%%%%%%%%%%%
        VIP_n=vVIP_n';
        VIP_n=VIP_n(ones(1,ncell_N),:);
        S_VIP_n=0+sum(VIP_n.*A,2)'.*sumal;
        beta1 = (0+S_VIP_n)'./(KD_n+(0+S_VIP_n))';


    %%% GABA calcs %%%%%%%%%%%%
        GABA=0.1+gGABA_n';
        GABA=GABA(ones(1,ncell_N),:);
        GABA_a=gGABA_a';
        GABA_a=GABA_a(ones(1,ncell_N),:);
        S_GABA=sum(GABA.*A1,2).*sumalGABA'; 
        S_GABAr=(GABA_a.*B1)';
        S_GABAp = 0+S_GABA - (S_GABA.*B1')  +S_GABAr;

        
    %%% Glu calcs %%%%%%%%%%%%   
        Glu=gGlu_a';
        Glu=Glu(ones(1,ncell_N),:);
        S_Glu_n=0+(Glu.*B2)';
        beta2= S_Glu_n'.^nGluR./(Km_GluR+S_Glu_n'.^nGluR);


    end



[vv]= FiringRates((Ca_n),S_GABAp, Fir,CC_n,BC_n,MP_n);


%%% About the intracelluar pathways

vo_n=(v_vo_n*BC_n.^n_vo_n./(K_vo_n+BC_n.^n_vo_n))+(vm_GluR.*beta2);
kk1_n=v_kk_n*CC_n.^n_kk_n./(K_kk_n+CC_n.^n_kk_n);
vv2_n =VM2_n *(Ca_n.^M2_n)./(K2_n^M2_n+Ca_n.^M2_n);
vv3_n =1*(VM3_n.*(Ca_store_n.^M3_n)./(KR_n^M3_n+Ca_store_n.^M3_n)).*(Ca_n.^pA_n)./(KA_n^pA_n+Ca_n.^pA_n);
vK_n = switz*(VMK_n.*(Ca_n)./(kMK_n+(Ca_n))+V_b_n*beta1./(k_b_n+beta1));
vsPc_n =vsP0_n+ CT_n.*CB_n./(KC_n+CB_n);


    for j=1:ncell_N

        dydt((j-1)*ns_N+1,1) = vo_n(j)+ 0.0003*IP3_n-vv2_n(j)+vv3_n(j)+kf_n*Ca_store_n(j) - kk1_n(j)*Ca_n(j).^(2); 
        dydt((j-1)*ns_N+2,1) = (vv2_n(j)-vv3_n(j)-kf_n*Ca_store_n(j));
        dydt((j-1)*ns_N+3,1) = vsPc_n(j)*BN_n(j)^n_n/(KAP_n^n_n+BN_n(j)^n_n)-vmP_n*MP_n(j)/(KmP_n+MP_n(j))-kdmp_n*MP_n(j);
        dydt((j-1)*ns_N+4,1) = vsC_n*BN_n(j)^n_n/(KAC_n^n_n+BN_n(j)^n_n)-vmC_n*MC_n(j)/(KmC_n+MC_n(j))-kdmc_n*MC_n(j);
        dydt((j-1)*ns_N+5,1) = vsB_n(j)*KIB_n^m_n/(KIB_n^m_n+BN_n(j)^m_n)-vmB_n(j)*MB_n(j)/(KmB_n+MB_n(j))-kdmb_n*MB_n(j);
        dydt((j-1)*ns_N+6,1) = ksP_n*MP_n(j)-V1P_n*PC_n(j)/(Kp_n+PC_n(j))+V2P_n*PCP_n(j)/(Kdp_n+PCP_n(j))+k4_n*PCC_n(j)-k3_n*PC_n(j)*CC_n(j)-kdn_n*PC_n(j);
        dydt((j-1)*ns_N+7,1) = ksC_n*MC_n(j)-V1C_n*CC_n(j)/(Kp_n+CC_n(j))+V2C_n*CCP_n(j)/(Kdp_n+CCP_n(j))+k4_n*PCC_n(j)-k3_n*PC_n(j)*CC_n(j)-kdnc_n*CC_n(j);
        dydt((j-1)*ns_N+8,1) = V1P_n*PC_n(j)/(Kp_n+PC_n(j))-V2P_n*PCP_n(j)/(Kdp_n+PCP_n(j))-vdPC_n*PCP_n(j)/(Kd_n+PCP_n(j))-kdn_n*PCP_n(j);
        dydt((j-1)*ns_N+9,1) = V1C_n*CC_n(j)/(Kp_n+CC_n(j))-V2C_n*CCP_n(j)/(Kdp_n+CCP_n(j))-vdCC_n*CCP_n(j)/(Kd_n+CCP_n(j))-kdn_n*CCP_n(j);
        dydt((j-1)*ns_N+10,1) = -V1PC_n*PCC_n(j)/(Kp_n+PCC_n(j))+V2PC_n*PCCP_n(j)/(Kdp_n+PCCP_n(j))-k4_n*PCC_n(j)+k3_n*PC_n(j)*CC_n(j)+k2_n*PCN_n(j)-k1_n*PCC_n(j)-kdn_n*PCC_n(j);
        dydt((j-1)*ns_N+11,1) = -V3PC_n*PCN_n(j)/(Kp_n+PCN_n(j))+V4PC_n*PCNP_n(j)/(Kdp_n+PCNP_n(j))-k2_n*PCN_n(j)+k1_n*PCC_n(j)-k7_n*BN_n(j)*PCN_n(j)+k8_n*IN_n(j)-kdn_n*PCN_n(j);
        dydt((j-1)*ns_N+12,1) = V1PC_n*PCC_n(j)/(Kp_n+PCC_n(j))-V2PC_n*PCCP_n(j)/(Kdp_n+PCCP_n(j))-vdPCC_n*PCCP_n(j)/(Kd_n+PCCP_n(j))-kdn_n*PCCP_n(j);
        dydt((j-1)*ns_N+13,1) = V3PC_n*PCN_n(j)/(Kp_n+PCN_n(j))-V4PC_n*PCNP_n(j)/(Kdp_n+PCNP_n(j))-vdPCN_n*PCNP_n(j)/(Kd_n+PCNP_n(j))-kdn_n*PCNP_n(j);
        dydt((j-1)*ns_N+14,1) = ksB_n*MB_n(j)-V1B_n*BC_n(j)/(Kp_n+BC_n(j))+V2B_n*BCP_n(j)/(Kdp_n+BCP_n(j))-k5_n*BC_n(j)+k6_n*BN_n(j)-kdn_n*BC_n(j);
        dydt((j-1)*ns_N+15,1) = V1B_n*BC_n(j)/(Kp_n+BC_n(j))-V2B_n*BCP_n(j)/(Kdp_n+BCP_n(j))-vdBC_n*BCP_n(j)/(Kd_n+BCP_n(j))-kdn_n*BCP_n(j);
        dydt((j-1)*ns_N+16,1) = -V3B_n*BN_n(j)/(Kp_n+BN_n(j))+V4B_n*BNP_n(j)/(Kdp_n+BNP_n(j))+k5_n*BC_n(j)-k6_n*BN_n(j)-k7_n*BN_n(j)*PCN_n(j)+k8_n*IN_n(j)-kdn_n*BN_n(j);
        dydt((j-1)*ns_N+17,1) = V3B_n*BN_n(j)/(Kp_n+BN_n(j))-V4B_n*BNP_n(j)/(Kdp_n+BNP_n(j))-vdBN_n*BNP_n(j)/(Kd_n+BNP_n(j))-kdn_n*BNP_n(j);
        dydt((j-1)*ns_N+18,1) = -k8_n*IN_n(j)+k7_n*BN_n(j)*PCN_n(j)-vdIN_n*IN_n(j)/(Kd_n+IN_n(j))-kdn_n*IN_n(j);
        dydt((j-1)*ns_N+19,1) = (vP_n/WT_n)*(vK_n(j)/vP_n*(1-CB_n(j))/(K_1_n+1-CB_n(j))-CB_n(j)/(K_2_n+CB_n(j)));
        dydt((j-1)*ns_N+20,1) = 0.5*vv(j).^1.9/(20+vv(j).^1.9)-0.5*vVIP_n(j).^0.2;
        dydt((j-1)*ns_N+21,1) = 0.5*vv(j).^1.9/(20+vv(j).^1.9)-0.5*gGABA_n(j).^0.2;

    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Astrocytes
      

%%Coupling


%%% no neurotransmitter release and coupling for the first 150h 

    if t<150
      switz=0;
   
      beta1_a=zeros(ncell_A,1);
      
    else 
        switz=1;
        
    %%% VIP calcs %%%%%%%%%%%%%

        VIP_a=vVIP_n';
        VIP_a=VIP_a(ones(1,ncell_A),:);
        S_VIP_a = (sum(VIP_a'.*B,1)).*sumal_VIP_NA;
        beta1_a = (0+S_VIP_a)'./(KD_a+(0+S_VIP_a))';

    end

%% 

%%% About the intracelluar pathways

vo_a=(v_vo_a*BC_a.^n_vo_a./(K_vo_a+BC_a.^n_vo_a));
kk1_a=v_kk_a*CC_a.^n_kk_a./(K_kk_a+CC_a.^n_kk_a);
vv2_a =VM2_a *(Ca_a.^M2_a)./(K2_a^M2_a+Ca_a.^M2_a);
vv3_a =1*(VM3_a.*(Ca_store_a.^M3_a)./(KR_a^M3_a+Ca_store_a.^M3_a)).*(Ca_a.^pA_a)./(KA_a^pA_a+Ca_a.^pA_a);
vK_a = switz*(VMK_a.*(Ca_a)./(kMK_a+(Ca_a))+V_b_a*beta1_a./(k_b_a+beta1_a));
vsPc_a =vsP0_a+ CT_a.*CB_a./(KC_a+CB_a);


%%% About GABA

v_GABA_a = vo_GABA_a + (vm_GABA_a*PCN_a.^n_GABA_a/(Km_GABA_a+PCN_a.^n_GABA_a));
v_iGAT = vm_iGAT.*CC_a.^n_GAT_a./(Km_GAT_a+CC_a.^n_GAT_a);
v_upGAT_a = (v_iGAT).*gGABA_a.^n_iGAT./(Km_iGAT+gGABA_a.^n_iGAT);
v_upGAT_n = (v_iGAT).*(S_GABA.*B1') .^n_iGAT./(Km_iGAT+(S_GABA.*B1') .^n_iGAT);
v_reGAT_a = (v_iGAT).*gGABA_store_a.^n_iGAT./(Km_iGAT+gGABA_store_a.^n_iGAT);

     
%%% About glutamate

v_Glu_a = vo_Glu_a + (vm_Glu_a.*CC_a.^n_Glu_a./(Km_Glu_a+CC_a.^n_Glu_a));
v_iEAAT = vm_iEAAT.*MC_a.^n_EAAT_a./(Km_EAAT_a+MC_a.^n_EAAT_a);
v_reEAAT_a =(v_iEAAT+0).*gGlu_store_a.^n_iEAAT./(Km_iEAAT+gGlu_store_a.^n_iEAAT);
v_upEAAT_a = (v_iEAAT+0).*gGlu_a.^n_iEAAT./(Km_iEAAT+gGlu_a.^n_iEAAT);

     
    for k=m
        

    dydt((ncell_N*ns_N)+(k-1)*ns_A+1,1) = vo_a(k)+ 0.0003*IP3_a-vv2_a(k)+vv3_a(k)+kf_a*Ca_store_a(k) - kk1_a(k)*Ca_a(k).^(2); 
    dydt((ncell_N*ns_N)+(k-1)*ns_A+2,1) = (vv2_a(k)-vv3_a(k)-kf_a*Ca_store_a(k)); 
    dydt((ncell_N*ns_N)+(k-1)*ns_A+3,1) = vsPc_a(k)*BN_a(k)^n_a/(KAP_a^n_a+BN_a(k)^n_a)-vmP_a*MP_a(k)/(KmP_a+MP_a(k))-kdmp_a*MP_a(k);                     
    dydt((ncell_N*ns_N)+(k-1)*ns_A+4,1) = vsC_a*BN_a(k)^n_a/(KAC_a^n_a+BN_a(k)^n_a)-vmC_a*MC_a(k)/(KmC_a+MC_a(k))-kdmc_a*MC_a(k);                   
    dydt((ncell_N*ns_N)+(k-1)*ns_A+5,1) = vsB_a*KIB_a^m_a/(KIB_a^m_a+BN_a(k)^m_a)-vmB_a*MB_a(k)/(KmB_a+MB_a(k))-kdmb_a*MB_a(k);                   
    dydt((ncell_N*ns_N)+(k-1)*ns_A+6,1) = ksP_a*MP_a(k)-V1P_a*PC_a(k)/(Kp_a+PC_a(k))+V2P_a*PCP_a(k)/(Kdp_a+PCP_a(k))+k4_a*PCC_a(k)-k3_a*PC_a(k)*CC_a(k)-kdn_a*PC_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+7,1) = ksC_a*MC_a(k)-V1C_a*CC_a(k)/(Kp_a+CC_a(k))+V2C_a*CCP_a(k)/(Kdp_a+CCP_a(k))+k4_a*PCC_a(k)-k3_a*PC_a(k)*CC_a(k)-kdnc_a*CC_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+8,1) = V1P_a*PC_a(k)/(Kp_a+PC_a(k))-V2P_a*PCP_a(k)/(Kdp_a+PCP_a(k))-vdPC_a*PCP_a(k)/(Kd_a+PCP_a(k))-kdn_a*PCP_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+9,1) = V1C_a*CC_a(k)/(Kp_a+CC_a(k))-V2C_a*CCP_a(k)/(Kdp_a+CCP_a(k))-vdCC_a*CCP_a(k)/(Kd_a+CCP_a(k))-kdn_a*CCP_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+10,1) = -V1PC_a*PCC_a(k)/(Kp_a+PCC_a(k))+V2PC_a*PCCP_a(k)/(Kdp_a+PCCP_a(k))-k4_a*PCC_a(k)+k3_a*PC_a(k)*CC_a(k)+k2_a*PCN_a(k)-k1_a*PCC_a(k)-kdn_a*PCC_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+11,1) = -V3PC_a*PCN_a(k)/(Kp_a+PCN_a(k))+V4PC_a*PCNP_a(k)/(Kdp_a+PCNP_a(k))-k2_a*PCN_a(k)+k1_a*PCC_a(k)-k7_a*BN_a(k)*PCN_a(k)+k8_a*IN_a(k)-kdn_a*PCN_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+12,1) = V1PC_a*PCC_a(k)/(Kp_a+PCC_a(k))-V2PC_a*PCCP_a(k)/(Kdp_a+PCCP_a(k))-vdPCC_a*PCCP_a(k)/(Kd_a+PCCP_a(k))-kdn_a*PCCP_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+13,1) = V3PC_a*PCN_a(k)/(Kp_a+PCN_a(k))-V4PC_a*PCNP_a(k)/(Kdp_a+PCNP_a(k))-vdPCN_a*PCNP_a(k)/(Kd_a+PCNP_a(k))-kdn_a*PCNP_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+14,1) = ksB_a*MB_a(k)-V1B_a*BC_a(k)/(Kp_a+BC_a(k))+V2B_a*BCP_a(k)/(Kdp_a+BCP_a(k))-k5_a*BC_a(k)+k6_a*BN_a(k)-kdn_a*BC_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+15,1) = V1B_a*BC_a(k)/(Kp_a+BC_a(k))-V2B_a*BCP_a(k)/(Kdp_a+BCP_a(k))-vdBC_a*BCP_a(k)/(Kd_a+BCP_a(k))-kdn_a*BCP_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+16,1) = -V3B_a*BN_a(k)/(Kp_a+BN_a(k))+V4B_a*BNP_a(k)/(Kdp_a+BNP_a(k))+k5_a*BC_a(k)-k6_a*BN_a(k)-k7_a*BN_a(k)*PCN_a(k)+k8_a*IN_a(k)-kdn_a*BN_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+17,1) = V3B_a*BN_a(k)/(Kp_a+BN_a(k))-V4B_a*BNP_a(k)/(Kdp_a+BNP_a(k))-vdBN_a*BNP_a(k)/(Kd_a+BNP_a(k))-kdn_a*BNP_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+18,1) = -k8_a*IN_a(k)+k7_a*BN_a(k)*PCN_a(k)-vdIN_a*IN_a(k)/(Kd_a+IN_a(k))-kdn_a*IN_a(k);
    dydt((ncell_N*ns_N)+(k-1)*ns_A+19,1) = (vP_a/WT_a)*(vK_a(k)/vP_a*(1-CB_a(k))/(K_1_a+1-CB_a(k))-CB_a(k)/(K_2_a+CB_a(k)));
    dydt((ncell_N*ns_N)+(k-1)*ns_A+20,1) =  v_reGAT_a(k)- v_upGAT_a(k)  - k_dGABAa*gGABA_a(k).^n_dGABAa ;
    dydt((ncell_N*ns_N)+(k-1)*ns_A+21,1) =  v_GABA_a(k) + v_upGAT_n(k)+ v_upGAT_a(k)- v_reGAT_a(k) - k_dGABAa_store*gGABA_a(k).^n_dGABAa_store ;
    dydt((ncell_N*ns_N)+(k-1)*ns_A+22,1) =  v_reEAAT_a(k) - v_upEAAT_a(k) - k_dGlua*gGlu_a(k).^n_dGlua;
    dydt((ncell_N*ns_N)+(k-1)*ns_A+23,1) =  v_Glu_a(k)+ v_upEAAT_a(k)- v_reEAAT_a(k) - k_dGlua_store*gGlu_store_a(k).^n_dGlua_store;


    end
    
    
end


