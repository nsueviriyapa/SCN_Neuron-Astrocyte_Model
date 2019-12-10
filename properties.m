function vv = properties (Ca_in,ENa,EK,q,k,T,Ca_out,Cl_out,Cl_in,F,BC)

Far =F(9);
Kout = F(22);
Naout = F(23);
PCa=F(29);
PCl= F(30);
PNa=F(32);
R=F(34);
vPK=F(46);
KPK=F(47);
nPK=F(48);

% Permeability 

PK=(vPK*BC.^nPK./(KPK+BC.^nPK)); 
thetaNa=exp(ENa*q/k/T/1000);
thetaK = exp(EK*q/k/T/1000);
Kin = Kout/thetaK;
Nain = Naout/thetaNa;
Cl_in=Cl_in';

% Spangler's equations

alpha = 4*PCa.*(Ca_in)*10^-3 + PK.*(Kin) + PNa.*(Nain)+PCl.*Cl_out;
bitaa = PK*(Kin)-PK*Kout + PNa*(Nain) - PNa*Naout +PCl*Cl_out-PCl*(Cl_in)';
c=-(4*PCa.*(Ca_out)*10^-3 + PK*(Kout) + PNa*(Naout)+PCl*(Cl_in)');
psi =  (-bitaa + sqrt(bitaa.^2 - 4*alpha.*c))./(2.*alpha);

% Membrane Potential

 vv= R*T/Far*2.303*log10(psi)*1000; 
 
 