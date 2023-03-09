clc 
clear
%% Body Partition's Dimension
% Main Body (Divided into three main parts: Fore, Middle, Aft
Mt = 52;
Lfb = .115;
Lmb = 1.245;
Lab = .24;
Ltotal = Lfb+Lmb+Lab;
D = .23;
R = D/2;
Vfb = 2/3*pi*(Lfb/2)^3;
Vmb = pi*(R)^2*Lmb;
Vab = 1/3*pi*(R)^2*Lab;

CoG = transpose([0 0 0]);
CoB = transpose([0 0 -.02]);

% Appendages
La1 = .023; % Cylinder
Da1 = .046;
Ra1 = Da1/2;
Va1 = pi*(Ra1)^2*La1;

La2 = .05; % Cylinder
Da2 = .046;
Ra2 = Da2/2;
Va2 = pi*(Ra2)^2*La2;

La3 = .069; % Square Prism
Ba3 = .0345;
Ha3 = .2415;
Va3 = La3*Ba3*Ha3;

La4 = .069; % Square Prism
Ba4 = .023;
Ha4 = .0115;
Va4 = La4*Ba4*Ha4;

% Thruster 
Lth = .231;
Dth = .1;
Rth = Dth/2;
Vth = pi*(Dth/2)^2*Lth;

%% Homogenous Density
Vb = Vfb+Vmb+Vab;
Va = Va1+Va2+Va3+Va4;

Vt = Vb+Va+(2*Vth);

rho = Mt/Vt;

%% Parts' Mass Calculation
% Main Body
Mfb = Vfb*rho;
Mmb = Vmb*rho;
Mab = Vab*rho;

% Appendages
Ma1 = Va1*rho;
Ma2 = Va2*rho;
Ma3 = Va3*rho;
Ma4 = Va4*rho;

% Thruster
Mth = Vth*rho;

%% Parts'Moment Inertia (In Respect to Its Own CoG)
% Main
Ifb = halfsphere_inertia(Mfb,R,Lfb);
Imb = cylinder_inertia(Mmb,R,Lmb);
Iab = cone_inertia(Mab, R, Lab); % In respect to the tip of the cone

% Appendages
Ia1 = cylinder_inertia(Ma1,Ra1,La1);
Ia2 = cylinder_inertia(Ma2,Ra2,La2);
Ia3 = squareprism_inertia(Ma3,La3,Ba3,Ha3);
Ia4 = squareprism_inertia(Ma4,La4,Ba4,Ha4);

% Thruster
Ith = cylinder_inertia(Mth,Rth,Lth);

%% Parts' Inertia Matrix (Its Own CoG to RB)
% Main Body
MIfbC = [mass_matrix(Mfb) zeros(3) ;
        zeros(3) Ifb];
MImbC = [mass_matrix(Mmb) zeros(3) ;
        zeros(3) Imb];
MIabC = [mass_matrix(Mab) zeros(3) ;
        zeros(3) Iab];

% Appendages
MIa1C = [zeros(3) zeros(3) ;
        zeros(3) zeros(3)];
MIa2C = [zeros(3) zeros(3) ;
        zeros(3) zeros(3)];
MIa3C = [zeros(3) zeros(3) ;
        zeros(3) zeros(3)];
MIa4C = [zeros(3) zeros(3) ;
        zeros(3) zeros(3)];

%Thruster
MIthC = [mass_matrix(Mth) zeros(3) ;
        zeros(3) Ith];

%% Parts' Inertia Matrix Transformation to RB CoG
% Part Position Defenition In Respect to RB [p,q,r]
Pfb = transpose([(3*R/8+.695) 0 0]);
Pmb = transpose([.0675 0 0]);
Pab = transpose([-(Lab/4+.56) 0 0]);
PabI = transpose([-.8 0 0]); % In respect to cone's tip

Pa1 = transpose([.638 0 La1/2+R]);
Pa2 = transpose([.371 0 La2/2+R]);
Pa3 = transpose([-.453 0 Ha3/2+R]);
Pa4 = transpose([-.555 0 Ha4/2+R]);

Pthl = transpose([-.59 -.17 0]);
Pthr = transpose([-.59 .17 0]);

% Inertia Matrix Transformation 
MIfb = HT_matrix(Pfb(1),Pfb(2),Pfb(3))*MIfbC*H_matrix(Pfb(1),Pfb(2),Pfb(3));
MImb = HT_matrix(Pmb(1),Pmb(2),Pmb(3))*MImbC*H_matrix(Pmb(1),Pmb(2),Pmb(3));
MIab = HT_matrix(PabI(1),PabI(2),PabI(3))*MIabC*H_matrix(PabI(1),PabI(2),PabI(3));

MIa1 = HT_matrix(Pa1(1),Pa1(2),Pa1(3))*MIa1C*H_matrix(Pa1(1),Pa1(2),Pa1(3));
MIa2 = HT_matrix(Pa2(1),Pa2(2),Pa2(3))*MIa2C*H_matrix(Pa2(1),Pa2(2),Pa2(3));
MIa3 = HT_matrix(Pa3(1),Pa3(2),Pa3(3))*MIa3C*H_matrix(Pa3(1),Pa3(2),Pa3(3));
MIa4 = HT_matrix(Pa4(1),Pa4(2),Pa4(3))*MIa4C*H_matrix(Pa4(1),Pa4(2),Pa4(3));

MIthl = HT_matrix(Pthl(1),Pthl(2),Pthl(3))*MIthC*H_matrix(Pthl(1),Pthl(2),Pthl(3));
MIthr = HT_matrix(Pthr(1),Pthr(2),Pthr(3))*MIthC*H_matrix(Pthr(1),Pthr(2),Pthr(3));

MI_Total = MIfb+MImb+MIab+MIthl+MIthr

%% Main Body Added Mass Matrix 
rhow = 1000;
% For m11 using spheroid method
a = (Lab+Lmb+Lfb)/2;
b = R;
e = sqrt(1-(b/a)^2);
alpha0 = 2*((1-e^2)/e^3)*((log((1+e)/(1-e)))/2-e);
beta0 = (1/e^2)-((1-e^2)/2*e^3)*log((1+e)/(1-e));
k1 = alpha0/(2-alpha0);
k2 = beta0/(2-beta0);
k4 = 0;
k5 = (e^4*(beta0-alpha0))/((2-e^2)*(2*e^2-(2-e^2)*(beta0-alpha0)));
m11 = 4/3*rho*pi*a*b^2*k1;

% Integration on x direction
CA2 = 1;
% Integration on the aft part
xl_ay = -.8;
xu_ay = -.56;
m22a = rhow*CA2*pi*0.479^2*((xu_ay^3/3 + 0.8*xu_ay^2 + 0.64*xu_ay)-(xl_ay^3/3 + 0.8*xl_ay^2 + 0.64*xl_ay));
m26a = rhow*CA2*pi*0.479^2*((xu_ay^4/4 + 0.53*xu_ay^3 + 0.32*xu_ay^2)-(xl_ay^4/4 + 0.53*xl_ay^3 + 0.32*xl_ay^2));
m66a = rhow*CA2*pi*0.479^2*((xu_ay^5/5 + 0.4*xu_ay^4 + 0.213*xu_ay^3)-(xl_ay^5/5 + 0.4*xl_ay^4 + 0.213*xl_ay^3)); 

% Integration on the middle part
xl_my = -.56;
xu_my = .685;
m22m = rhow*CA2*pi*0.013225*(xu_my - xl_my);
m26m = rhow*CA2*pi*0.0066125*(xu_my^2 - xl_my^2);
m66m = rhow*CA2*pi*0.004408*(xu_my^3 - xl_my^3);

% Integration on the fore part
xl_fy = .685;
xu_fy = .8;
m22f = rhow*CA2*pi*((-xu_fy^3/3 + .685*xu_fy^2 - .456*xu_fy) - (-xl_fy^3/3 + .685*xl_fy^2 - .456*xl_fy));
m26f = rhow*CA2*pi*((-xu_fy^4/4 + .457*xu_fy^3 - .228*xu_fy^2) - (-xl_fy^4/4 + .457*xl_fy^3 - .228*xl_fy^2));
m66f = rhow*CA2*pi*((-xu_fy^5/5 + .3425*xu_fy^4 - .152*xu_fy^3) - (-xl_fy^5/5 + .3425*xl_fy^4 - .152*xl_fy^3));

% Integration on z direction
CA3 = 1;
% Integration of the aft part 1
xl_az1 = -.8;
xu_az1 = -.7055;
m33a1 = rhow*CA3*pi*0.479^2*((xu_az1^3/3 + 0.8*xu_az1^2 + 0.64*xu_az1) - (xl_az1^3/3 + 0.8*xl_az1^2 + 0.64*xl_az1));
m35a1 = -rhow*CA3*pi*0.479^2*((xu_az1^4/4 + 0.53*xu_az1^3 + 0.32*xu_az1^2) - (xl_az1^4/4 + 0.53*xl_az1^3 + 0.32*xl_az1^2));
m55a1= rhow*CA3*pi*0.479^2*((xu_az1^5/5 + 0.4*xu_az1^4 + 0.213*xu_az1^3) - (xl_az1^5/5 + 0.4*xl_az1^4 + 0.213*xl_az1^3));

% Integration of the aft part 2
Dmin = .08269; % Diameter at x = -.7055 using Rhinoceros 3D
Rmin = Dmin/2;
spanfin = Dth+R;

xl_az2 = -.7055;
xu_az2 = -.4745;

CA32 = 1-(R/spanfin)^2+(R/spanfin)^4;
m33a2 = rhow*CA3*pi*spanfin^2*(xu_az2 - xl_az2);
m35a2 = -rhow*CA3*pi*spanfin^2*(xu_az2^2/2 - xl_az2^2/2);
m55a2 = rhow*CA3*pi*spanfin^2*(xu_az2^3/3 - xl_az2^3/3);

alpha = asin((2*R*spanfin)/(spanfin^2+R^2));
f_alpha = 2*alpha^2-alpha*sin(4*alpha) + 0.5*sin(2*alpha);
m44a2 = rhow*R^4*((csc(alpha))^4*f_alpha-pi^2)/(2*pi);

% Integration on the middle part
xl_mz = -.4745;
xu_mz = .685;
m33m = rhow*CA3*pi*0.013225*(xu_mz - xl_mz);
m35m = -rhow*CA3*pi*0.0066125*(xu_mz^2 - xl_mz^2);
m55m = rhow*CA3*pi*0.004408*(xu_mz^3 - xl_mz^3);

% Integration on the front part
xl_fz = .685;
xu_fz = .8;
m33f = rhow*CA2*pi*((-xu_fz^3/3 + .685*xu_fz^2 - .456*xu_fz) - (-xl_fz^3/3 + .685*xl_fz^2 - .456*xl_fz));
m35f = -rhow*CA2*pi*((-xu_fz^4/4 + .457*xu_fz^3 - .228*xu_fz^2) - (-xl_fz^4/4 + .457*xl_fz^3 - .228*xl_fz^2));
m55f = rhow*CA2*pi*((-xu_fz^5/5 + .3425*xu_fz^4 - .152*xu_fz^3) - (-xl_fz^5/5 + .3425*xl_fz^4 - .152*xl_fz^3));

% Element list
m22 = m22a + m22m + m22f;
m26 = m26a + m26m + m26f;
m33 = m33a1 + m33a2 + m33m + m33f;
m35 = m35a1 + m35a2 + m35m + m35f;
m44 = m44a2;
m55 = m55a1 + m55a2 + m55m + m55f;
m66 = m66a + m66m + m66f;

% Added Mass Matrix Transformation to CoG
AMM = [m11 0 0 0 0 0 ;
        0 m22 0 0 0 m26;
        0 0 m33 0 m35 0;
        0 0 0 m44 0 0;
        0 0 m35 0 m55 0;
        0 m26 0 0 0 m66];

%% Appendages Added Mass Matrix
%AMM1 (Small Cylinder)
ma11_1 = rhow*1/2*4/3*pi*Ra1^3; % Using Added mass for 3D Object method (spheres)
ma22_1 = ma11_1;
ma33_1 = 0; % Object is not slender enough
maI_1 = spheroid_AMM(La1, Ra1);
AMM1BT = [ma11_1 0 0 0 0 0;
        0 ma22_1 0 0 0 0;
        0 0 ma33_1 0 0 0;
        0 0 0 maI_1(4) 0 0;
        0 0 0 0 maI_1(5) 0;
        0 0 0 0 0 maI_1(6)];

%AMM2 (Small Cylinder)
ma11_2 = rhow*1/2*4/3*pi*Ra2^3; % Using Added mass for 3D Object method (spheres)
ma22_2 = ma11_2;
ma33_2 = 0; % Object is not slender enough
maI_2 = spheroid_AMM(La2, Ra2);
AMM2BT = [ma11_2 0 0 0 0 0;
        0 ma22_2 0 0 0 0;
        0 0 ma33_2 0 0 0;
        0 0 0 maI_2(4) 0 0;
        0 0 0 0 maI_2(5) 0;
        0 0 0 0 0 maI_2(6)];

% AMM3 (Long Square Prism)
CARP311 = .934; % CA for b/a = 7, bo/ao = 7
CARP322 = .872; % CA for b/a = 4, bo/ao = 3.5
CASP3 = .19; % CA for b/a = 4, bo/ao = 3.5
ma11_3 = rhow*CARP311*pi*Ba3^2*Ha3; % Using Added mass for 3D Object method (flat rectangular plates)
ma22_3 = rhow*CARP322*pi*La3^2*Ha3; % Using Added mass for 3D Object method (flat rectangular plates)
ma33_3 = rhow*CASP3*La3^2*Ha3; % Object is not slender enough
maI_3 = spheroid_AMM(Ha3, La3);
AMM3BT = [ma11_3 0 0 0 0 0;
        0 ma22_3 0 0 0 0;
        0 0 ma33_3 0 0 0;
        0 0 0 maI_3(4) 0 0;
        0 0 0 0 maI_3(5) 0;
        0 0 0 0 0 maI_3(6)];

% AMM4 (Long Square Prism)
CARP411 = .579; % CA for b/a = 1,  bo/ao < 1
CARP422 = .579; % CA for b/a = 1  bo/ao <1
CASP4 = .68; % CA for b/a = 1 bo/ao <1
ma11_4 = rhow*CARP411*pi*Ba4^2*Ha4; % Using Added mass for 3D Object method (flat rectangular plates)
ma22_4 = rhow*CARP422*pi*La4^2*Ha4; % Using Added mass for 3D Object method (flat rectangular plates)
ma33_4 = rhow*CASP4*La4^2*Ha4; % Object is not slender enough
maI_4 = spheroid_AMM(Ha4, La4);
AMM4BT = [ma11_4 0 0 0 0 0;
        0 ma22_4 0 0 0 0;
        0 0 ma33_4 0 0 0;
        0 0 0 maI_4(4) 0 0;
        0 0 0 0 maI_4(5) 0;
        0 0 0 0 0 maI_4(6)];
%% Parts' Added Mass Matrix Transformation to RB CoG
% Transformation for Main Body Added Mass Matrix to RB CoG
PAM = transpose([0 0 -.02]); % Position of CoB of Main Body
AMMG = HT_matrix(PAM(1),PAM(2),PAM(3))*AMM*H_matrix(PAM(1),PAM(2),PAM(3));

% Transformation for Appendages Added Mass Matrix to RB CoG using each
% part's CoG
AMM1 = HT_matrix(Pa1(1),Pa1(2),Pa1(3))*AMM1BT*H_matrix(Pa1(1),Pa1(2),Pa1(3));
AMM2 = HT_matrix(Pa2(1),Pa2(2),Pa2(3))*AMM2BT*H_matrix(Pa2(1),Pa2(2),Pa2(3));
AMM3 = HT_matrix(Pa3(1),Pa3(2),Pa3(3))*AMM3BT*H_matrix(Pa3(1),Pa3(2),Pa3(3));
AMM4 = HT_matrix(Pa4(1),Pa4(2),Pa4(3))*AMM4BT*H_matrix(Pa4(1),Pa4(2),Pa4(3));

% Total Added Mass Matrix
AMM_Total = AMMG 

%% Velocity Transform
H_fb = H_matrix(Pfb(1),Pfb(2),Pfb(3));
H_mb = H_matrix(Pmb(1),Pmb(2),Pmb(3));
H_ab = H_matrix(Pab(1),Pab(2),Pab(3));

H_Mb = H_matrix(CoB(1),CoB(2),CoB(3));

H_a1 = H_matrix(Pa1(1),Pa1(2),Pa1(3));
H_a2 = H_matrix(Pa2(1),Pa2(2),Pa2(3));
H_a3 = H_matrix(Pa3(1),Pa3(2),Pa3(3));
H_a4 = H_matrix(Pa4(1),Pa4(2),Pa4(3));

H_thl = H_matrix(Pthl(1),Pthl(2),Pthl(3));
H_thr = H_matrix(Pthr(1),Pthr(2),Pthr(3));

% Part' Velocities
% velo_DVL = transpose([1 1 1 0 0 0]);

% velo_fb = H_fb.*velo_DVL;
% velo_mb = H_mb.*velo_DVL;
% Velo_ab = H_ab.*velo_DVL;
% 
% velo_Mb = H_Mb.*velo_DVL;
% 
% velo_a1 = H_a1.*velo_DVL;
% velo_a2 = H_a2.*velo_DVL;
% velo_a3 = H_a3.*velo_DVL;
% velo_a4 = H_a4.*velo_DVL;
% 
% velo_thl = H_thl.*velo_DVL;
% velo_thr = H_thr.*velo_DVL;

%% Main Body Drag Matrix
% Main Body x-direction
Cd11Mb = .4;
SxMb = pi*(D^2)/4;

% Main Body y-direction
Cd22Mb = 1.2;
Afb = pi*(R^2)/2;
Amb = Lmb*D;
Aab = Lab*D/2; 
SyMb = Afb+Amb+Aab;

% Main Body z-direction
Cd33Mb = 1.2;
SzMb = SyMb;

% Drag Matrix
HT_Mb = HT_matrix(CoB(1),CoB(2),CoB(3));
K_Mb = K_matrix(Ltotal, SxMb, SyMb, SzMb, D, D, Cd11Mb, Cd22Mb, Cd33Mb)
% Tau_Mb = HT_Mb*K_Mb*abs(velo_Mb)*(velo_Mb);

%% Appendage 1 Drag Matrix
% Main Body x-direction
Cd11a1 = .55; % with La1/Da1 = 0.5 (Finite Cylinder Vertical, extrapolation)
Sxa1 = La1*Da1;

% Main Body y-direction
Cd22a1 = .55; % with La1/Da1 = 0.5 (Finite Cylinder Vertical, extrapolation)
Sya1 = La1*Da1;

% Main Body z-direction
Cd33a1 = 1.1; % with La1/Da1 = 0.5 (Finite Cylinder Horizontal)
Sza1 = pi*(Da1^2)/4;

% Drag Matrix
HT_a1 = HT_matrix(Pa1(1),Pa1(2),Pa1(3));
K_a1 = K_matrix(La1, Sxa1, Sya1, Sza1, Da1, Da1, Cd11a1, Cd22a1, Cd33a1)
% Tau_a1 = HT_a1*K_a1*abs(velo_a1)*(velo_a1);

%% Appendage 2 Drag Matrix
% Main Body x-direction
Cd11a2 = 0.6; % with La2/Da2 ~ 1 (Finite Cylinder Vertical)
Sxa2 = La2*Da2;

% Main Body y-direction
Cd22a2 = 0.6; % with La2/Da2 ~ 1 (Finite Cylinder Vertical)
Sya2 = La2*Da2;

% Main Body z-direction
Cd33a2 = 0.9; % with La2/Da2 ~ 1 (Finite Cylinder Horizontal)
Sza2 = pi*(Da2^2)/4;

% Drag Matrix
HT_a2 = HT_matrix(Pa2(1),Pa2(2),Pa2(3));
K_a2 = K_matrix(La2, Sxa2, Sya2, Sza2, Da2, Da2, Cd11a2, Cd22a2, Cd33a2)
% Tau_a2 = HT_a2*K_a2*abs(velo_a2)*(velo_a2);

%% Appendage 3 Drag Matrix
% Main Body x-direction
Cd11a3 = 2.2; % with La3/Ha3 ~ 0.3 (Rectangular Rod, Interpolate)
Sxa3 = Ba3*Ha3;

% Main Body y-direction
Cd22a3 = 1.9; % with Ba3/Ha3 ~ 0.1 (Rectangular Rod)
Sya3 = La3*Ha3;

% Main Body z-direction
Cd33a3 = 1.3; % with Ha3/La3 ~ 3 (Rectangular Rod, Interpolate) 
Sza3 = Ba3*La3;

% Drag Matrix
HT_a3 = HT_matrix(Pa3(1),Pa3(2),Pa3(3));
K_a3 = K_matrix(Ha3, Sxa3, Sya3, Sza3, La3, Ba3, Cd11a3, Cd22a3, Cd33a3)
% Tau_a3 = HT_a3*K_a3*abs(velo_a3)*(velo_a3);

%% Appendage 4 Drag Matrix
% Main Body x-direction
Cd11a4 = 1.3; % with La4/Ha4 ~ 3 (Rectangular Rod, Interpolatiom)
Sxa4 = Ba4*Ha4;

% Main Body y-direction
Cd22a4 = 1.7; % with Ba4/Ha4 ~ 2 (Rectangular Rod)
Sya4 = La4*Ha4;

% Main Body z-direction
Cd33a4 = 1.9; % with Ha4/La4 ~ 0.1 (Rectangular Rod, Interpolation) 
Sza4 = Ba4*La4;

% Drag Matrix
HT_a4 = HT_matrix(Pa4(1),Pa4(2),Pa4(3));
K_a4 = K_matrix(Ha4, Sxa4, Sya4, Sza4, La4, Ba4, Cd11a4, Cd22a4, Cd33a4)
% Tau_a4 = HT_a4*K_a4*abs(velo_a4)*(velo_a4);

%% Left Thruster Drag Matrix
% Main Body x-direction
Cd11thl = 1.1; % with Lth/Dth ~ 0.5 (Finite Cylinder Horizontal, Extrapolation)
Sxthl = pi*(Dth^2)/4;

% Main Body y-direction
Cd22thl = 0.6; % with Lth/Dth ~ 1 (Finite Cylinder Vertical, Extrapolation)
Sythl = Lth*Dth;

% Main Body z-direction
Cd33thl = 0.6; % with Lth/Dth ~ 1 (Finite Cylinder Vertical, Extrapolation)
Szthl = Lth*Dth;

% Drag Matrix
HT_thl = HT_matrix(Pthl(1),Pthl(2),Pthl(3));
K_thl = K_matrix(Lth, Sxthl, Sythl, Szthl, Dth, Dth, Cd11thl, Cd22thl, Cd33thl)
% Tau_thl = HT_thl*K_thl*abs(velo_thl)*(velo_thl);

%% Right Thruster Drag Matrix
% Main Body x-direction
Cd11thr = 1.1; % with Lth/Dth ~ 0.5 (Finite Cylinder Horizontal, Extrapolation)
Sxthr = pi*(Dth^2)/4;

% Main Body y-direction
Cd22thr = 0.6; % with Lth/Dth ~ 1 (Finite Cylinder Vertical, Extrapolation)
Sythr = Lth*Dth;

% Main Body z-direction
Cd33thr = 0.6; % with Lth/Dth ~ 1 (Finite Cylinder Vertical, Extrapolation)
Szthr = Lth*Dth;

% Drag Matrix
HT_thr = HT_matrix(Pthr(1),Pthr(2),Pthr(3));
K_thr = K_matrix(Lth, Sxthr, Sythr, Szthr, Dth, Dth, Cd11thr, Cd22thr, Cd33thr)
% Tau_thr = HT_thr*K_thr*abs(velo_thr)*(velo_thr);

%% Total Drag Matrix
% Drag_Total = Tau_Mb + Tau_a1 + Tau_a2 + Tau_a3 + Tau_a4 + Tau_thl +Tau_thr;

filename = 'SparusWA.mat';
save(filename)



