function [AccG] = RovModel(Thrust,PosE,VitB)

global Para

filename1 = 'Sparus.mat';
filename2 = 'SparusWA.mat';
load(filename2)

%% Attitudes in earth frame
% z=PosE(3,1);
phi     = PosE(4,1)	;
theta   = PosE(5,1)	;

%% Gravity Force

Fg = 1* [-Para.P * sin(theta) ;
        Para.P * cos(theta)*sin(phi) ;
        Para.P * cos(theta)*cos(phi) ;
        0 ;
        0 ;
        0 ];
    
% Expressed in b and computed in G
    
%% Force d'Archimède

Fa_F = [Para.B * sin(theta) ;
        -Para.B * cos(theta)*sin(phi) ;
        -Para.B * cos(theta)*cos(phi) ;
        ];
%  Expressed in b


Fa_M = S_(Para.rb-Para.rg) * Fa_F ; % Computed in G

Fa = [ Fa_F ; Fa_M ] ;
%  Expressed in b and computed in G

%% Force de Coriolis

u = VitB(1,1)   ;
v = VitB(2,1)   ;
w = VitB(3,1)   ;
p = VitB(4,1)   ;   %Body fixed velocity roll (rad*s^(-1))
q = VitB(5,1)   ;   %Body fixed velocity pitch (rad*s^(-1))
r = VitB(6,1)   ;   %Body fixed velocity yaw (rad*s^(-1))

V_ = [u;v;w]    ;
W_ = [p;q;r]    ;  %General vector


% Wb :
Wb = [  S_(W_)       zeros(3,3) ;
        zeros(3,3)      S_(W_)       ];

M_11 = Para.Mg(1:3,1:3);
M_12 = Para.Mg(1:3,4:6);
M_21 = Para.Mg(4:6,1:3);
M_22 = Para.Mg(4:6,4:6);

% General coriolis matrix :
C_all = [ zeros(3,3)                          -S_(M_11*V_+M_12*W_) ;
          -S_(M_11*V_+M_12*W_)      -S_(M_21*V_+M_22*W_)];

%coriolis Force :
Fc = C_all * VitB   ;

%% Friction forces
% Veloctiy Transformation
Vit_0=VitB;
velo_Mb = H_Mb*Vit_0;
Ff_0 =  HT_Mb*Para.S0.Kq * abs(velo_Mb).*velo_Mb ;

Vit_1=VitB;
velo_a1 = H_a1*Vit_1;
Ff_1 =  HT_a1*Para.S1.Kq * abs(velo_a1).*velo_a1 ;

Vit_2=VitB;
velo_a2 = H_a2*Vit_2;
Ff_2 =  HT_a2*Para.S2.Kq * abs(velo_a2).*velo_a2 ; 

Vit_3=VitB;
velo_a3 = H_a3*Vit_3;
Ff_3 =  HT_a3*Para.S3.Kq * abs(velo_a3).*velo_a3 ;

Vit_4=VitB;
velo_a4 = H_a4*Vit_4;
Ff_4 =  HT_a4*Para.S4.Kq * abs(velo_a4).*velo_a4 ;

Vit_tl=VitB;
velo_tl = H_thl*Vit_tl;
Ff_tl =  HT_thl*Para.Stl.Kq * abs(velo_tl).*velo_tl ;

Vit_tr=VitB;
velo_tr = H_thr*Vit_tr;
Ff_tr =  HT_thr*Para.Str.Kq * abs(velo_tr).*velo_tr ;

%% Propulsions Forces
Fp = Para.Eb * Thrust ;


%% Accelearion computation :
AccG = Para.Mg\(Ff_0+Ff_tl+Ff_1+Ff_2 +Ff_3+Ff_4+Ff_tr+Fa + Fg+ Fp- Fc) ; % Mg\ = Mg^-1 computed at the gravity center of the Sparus

