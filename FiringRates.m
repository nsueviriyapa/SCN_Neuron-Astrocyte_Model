function [v]= FiringRates(Ca_in,GABA,F,CC,BC,MP) 

% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ca_out = F(1);
Cl_ex=F(2);
Cl_o=F(3);
Cp=F(4);
E_ex=F(5);
EK_o=F(6);
EL_o=F(7);
ENa_o=F(8);%
g_inhib= F(10);
gKo=F(11);
gNa =F(12);
k = F(13);
K_R=F(14);
KCa= F(15);
KCl1=F(16);
KCl2=F(17);
Kex1=F(18);
Kex2=F(19);
Kgk=F(20);
KKCa=F(21);
nca=F(24);
nCl=F(25);
nex1=F(26);
nex2=F(27);
nKCa=F(28);
q = F(33);
T = F(35);
T_room= F(36);
V_R=F(37);
V_theta=F(38);
vCa=F(39);
vCl1=F(40);
vCl2=F(41);
Vex1= 1*F(42);
Vex2=1*F(43);
vgk=F(44);
vKCa=F(45);


% REVERSAL POTENTIALS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENa = ENa_o * T/(T_room);%
EK = (EK_o)*T/(T_room); %
EL =EL_o* T/(T_room); 
ECa =k*T/(2*q)*log(Ca_out./Ca_in)*1000; %

Cl_in=Cl_o+(MP./(KCl1+MP)*vCl1)+(GABA.^nCl./(KCl2+GABA.^nCl))*vCl2*1;
E_inhib = -k*T/(q)*log(Cl_ex./Cl_in)*1000';


% MEMBRANE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vrest =properties (Ca_in,ENa,EK,q,k,T,Ca_out,Cl_ex,Cl_in,F,BC); % 
theta =Vrest + V_theta ; 
Vreset= Vrest+4;
R=V_R*(Vrest)./(K_R + Vrest);  


% CONDUCTANCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INa = gNa.*(Vrest-ENa);
gK=(gKo+MP./(Kgk+MP)*vgk);
g_ex=(Vex1*abs(INa).^nex1./(Kex1+abs(INa).^nex1)+ ((Ca_in).^nex2)./(Kex2+(Ca_in).^nex2).*Vex2);
gL = 1./R ;
gCa=vCa*(MP.^nca./(KCa+MP.^nca));
gKCa =vKCa*(CC.^nKCa./(KKCa+CC.^nKCa));
I_inhib = g_inhib.*(Vrest- E_inhib);

%FINAL CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_star = (-g_inhib.*E_inhib -g_ex.*E_ex+gNa.*ENa + gCa.*ECa + gK.*EK + gL.*EL+ gKCa.*EK) ;
R_star = 1./(gNa+ gK + gL + gCa + gKCa- g_inhib - g_ex) ;
tau_m =Cp.*R_star;
v= -(tau_m.*log( (theta-(R_star.*I_star))./(Vreset-R_star.*I_star))).^-1 ;
