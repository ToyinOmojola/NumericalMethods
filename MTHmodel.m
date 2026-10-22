function ydot = MTHmodel(Tau, mf)
% Tau is the contact time, and mf is the mass fraction

% This function file was developed by Toyin Omojola in October 2026
% based on his thesis work at the University of Bath that was submitted in June 2019.
% The code is for steady-state fixed bed kinetic modelling
% Methanol-to-hydrocarbon conversion is used as an example but it can be
% adapted to examples in heterogeneous catalysis



%% globalise parameters
global k1fValue 

%% constants
%presssure (bar)
P_Tot = 1;

% total number of moles
% mf(i) is the mass fraction of specie of specie i
% current code has 33 species 

N_tot = (mf(1)/32)+(mf(2)/46)+(mf(3)/18)+(mf(4)/16)+(mf(5)/2)+(mf(6)/28)+(mf(7)/30)...
+(mf(8)/28)+(mf(9)/42)+(mf(10)/56)+(mf(11)/70)+(mf(12)/84)+(mf(13)/98)+(mf(14)/112)...
+(mf(15)/126)+(mf(16)/140)+(mf(17)/30)+(mf(18)/44)+(mf(19)/58)+(mf(20)/72)+(mf(21)/86)...
+(mf(22)/100)+(mf(23)/114)+(mf(24)/128)+(mf(25)/142)+(mf(26)/76)+(mf(27)/90)...
+(mf(28)/78)+(mf(29)/92)+(mf(30)/106)+(mf(31)/120)+(mf(32)/134)+(mf(33)/28);

%% Define Adsorption Equilibrium Constants 
% Ki is the adsorption equilibrium constant of specie i
% adsorption constant can be allocated according to functional groups

K1 = 190;
K2 = 250;
K3 = 209;
K4 = 0.005;
K5 = 0.01;
K6 = (1/3)*K4;
K7 = 400;
K8 = 1.5;
K9 = 1.3;
K10 = 0.7;
K11 = 0.3;
K12 = 1.00;
K13 = 1.53;
K14 = 1.53;
K15 = 1.53;
K16 = 1.53;
K17 = K4;
K18 = K4;
K19 = K4;
K20 = K4;
K21 = K4;
K22 = K4;
K23 = K4;
K24 = K4;
K25 = K4;
K26 = K2; 
K27 = K2;
K28=7.81;
K29=7.81;
K30=4;
K31=7.81;
K32=7.81;

%% Mole fraction definition
% specie (i) is defined 
x_CH3OH = (mf(1)/32)/N_tot;
x_CH3OCH3 = (mf(2)/46)/N_tot;
x_H2O = (mf(3)/18)/N_tot;
x_CH4=(mf(4)/16)/N_tot;
x_H2=(mf(5)/2)/N_tot;
x_CO=(mf(6)/28)/N_tot;
x_CH2O=(mf(7)/30)/N_tot;
x_C2H4 = (mf(8)/28)/N_tot;
x_C3H6 = (mf(9)/42)/N_tot;
x_C4H8 = (mf(10)/56)/N_tot;
x_C5H10 = (mf(11)/70)/N_tot;
x_C6H12 = (mf(12)/84)/N_tot;
x_C7H14 = (mf(13)/98)/N_tot;
x_C8H16 = (mf(14)/112)/N_tot;
x_C9H18 = (mf(15)/126)/N_tot;
x_C10H20 = (mf(16)/140)/N_tot;
x_C2H6 = (mf(17)/30)/N_tot;
x_C3H8 = (mf(18)/44)/N_tot;
x_C4H10 = (mf(19)/58)/N_tot;
x_C5H12 = (mf(20)/72)/N_tot;
x_C6H14 = (mf(21)/86)/N_tot;
x_C7H16 = (mf(22)/100)/N_tot;
x_C8H18 = (mf(23)/114)/N_tot;
x_C9H20 = (mf(24)/128)/N_tot;
x_C10H22 = (mf(25)/142)/N_tot;
x_C3H8O2 = (mf(26)/76)/N_tot;
x_C4H10O2 = (mf(27)/90)/N_tot;
x_C6H6 = (mf(28)/78)/N_tot;
x_C7H8 = (mf(29)/92)/N_tot;
x_C8H10 = (mf(30)/106)/N_tot;
x_C9H12 = (mf(31)/120)/N_tot;
x_C10H14 = (mf(32)/134)/N_tot;
x_N2 = (mf(33)/28)/N_tot;

%% Define fraction of vacant sites
% Fraction of vacant sites is based on adsorption of all species;
% assumptions can be made

HZ=1/(1+(K1*x_CH3OH*P_Tot)+(K2*x_CH3OCH3*P_Tot)+(K3*x_H2O*P_Tot)+(K4*x_CH4*P_Tot)+...
    (K5*x_H2*P_Tot)+(K6*x_CO*P_Tot)+(K7*x_CH2O*P_Tot)+(K8*x_C2H4*P_Tot)+(K9*x_C3H6*P_Tot)+...
    (K10*x_C4H8*P_Tot)+(K11*x_C5H10*P_Tot)+(K12*x_C6H12*P_Tot)+(K13*x_C7H14*P_Tot)+...
    (K14*x_C8H16*P_Tot)+(K15*x_C9H18*P_Tot)+(K16*x_C10H20*P_Tot)+(K17*x_C2H6*P_Tot)+...
    (K18*x_C3H8*P_Tot)+(K19*x_C4H10*P_Tot)+(K20*x_C5H12*P_Tot)+(K21*x_C6H14*P_Tot)+...
    (K22*x_C7H16*P_Tot)+(K23*x_C8H18*P_Tot)+(K24*x_C9H20*P_Tot)+(K25*x_C10H22*P_Tot)+...
(K26*x_C3H8O2*P_Tot)+(K27*x_C4H10O2*P_Tot)+(K28*x_C6H6*P_Tot)+(K29*x_C7H8*P_Tot)+...
(K30*x_C8H10*P_Tot)+(K31*x_C9H12*P_Tot)+(K32*x_C10H14*P_Tot)+(K33*x_N2*P_Tot));

%% reallocate globalised parameters 
k1f=K1fValue;

%% define rate equations
% rate equations
r1f=k1f*K1*(x_CH3OH^2)*(P_Tot^2)*HZ;
r1b=k1b*K2*x_CH3OCH3*x_H2O*(P_Tot^2)*HZ;

%% Differential equations for 33 species
ydot=zeros(33,1);

ydot(1)=32*((-2*r1f)+(2*r1b)-r2+rox2-rox3-r4d-r8-r9-r10-r11-r12-r13-r14-r15+r16+r17+r18+r19+r20+r21+r22+r23-r40-r41-r42-r43+r44+r45+r46+r47);
ydot(2)=46*(r1f-r1b-r3-(2*rox1)-rox2-r4d-r5d-r16-r17-r18-r19-r20-r21-r22-r23-r44-r45-r46-r47);
ydot(3)=18*(r1fr1b+rox3+r4id+(2*r4d)+r5d+r8+r9+r10+r11+r12+r13+r14+r15+r40+r41+r42+r43);

end 
