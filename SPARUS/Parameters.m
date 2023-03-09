function Para=Parameters()
global Para

filename1 = 'Sparus.mat';
filename2 = 'SparusWA.mat';
load(filename2)

%% Initial Speed and position in Earth-fixed frame

Para.ICPos = [0 0 2 0 0 0];
Para.ICSpeed = [0 0 0 0 0 0] ;

%% General parameters
Para.rho_water = 1000 ;                     % Masse volumique de l'eau (kg/m^3)
Para.R = 0.115 ;                             % Sparus Radius (m)
Para.D = Para.R * 2;
Para.L = 1.6;  	                            % Sparus length (m)
Para.Lfb = 0.115;
Para.Lmb = 1.245;
Para.Lab = 0.24;
Para.Ltotal = Para.Lfb+Para.Lmb+Para.Lab;
Para.Vfb = 2/3*pi*(Para.Lfb/2)^3;
Para.Vmb = pi*(Para.R)^2*Para.Lmb;
Para.Vab = 1/3*pi*(Para.R)^2*Para.Lab;
Para.m = 52 ; 	                            % Sparus mass (kg)
Para.mb = 52.1;                           	% Sparus buoyancy mass  (kg) 
Para.g = 9.81 ;                             % Earth Gravity (m*s^(-2))
Para.P = Para.m * Para.g;	                % Sparus weight (N)
Para.B = Para.mb * Para.g;	                % Buoyancy (N)

% Appendages
Para.La1 = .023; % Cylinder
Para.Da1 = .046;
Para.Ra1 = Para.Da1/2;
Para.Va1 = pi*(Para.Ra1)^2*Para.La1;

Para.La2 = .05; % Cylinder
Para.Da2 = .046;
Para.Ra2 = Para.Da2/2;
Para.Va2 = pi*(Para.Ra2)^2*Para.La2;

Para.La3 = .069; % Square Prism
Para.Ba3 = .0345;
Para.Ha3 = .2415;
Para.Va3 = Para.La3*Para.Ba3*Para.Ha3;

Para.La4 = .069; % Square Prism
Para.Ba4 = .023;
Para.Ha4 = .0115;
Para.Va4 = Para.La4*Para.Ba4*Para.Ha4;

% Thruster 
Para.Lth = .0231;
Para.Dth = .1;
Para.Rth = Para.Dth/2;
Para.Vth = pi*(Para.Dth/2)^2*Para.Lth;

% Rho Calculation
Para.Vb = Para.Vfb+Para.Vmb+Para.Vab;
Para.Va = Para.Va1+Para.Va2+Para.Va3+Para.Va4;

Para.Vt = Para.Vb+Para.Va+(2*Para.Vth);

Para.rho = Para.m/Para.Vt;

%% Center of gravity and Buoyancy position in body-fied frame

Para.xg = 0 ;    %x-positon of Center of gravity
Para.yg = 0 ;    %y-positon of Center of gravity
Para.zg = 0 ;    %z-positon of Center of gravity

Para.rg = [Para.xg Para.yg Para.zg]' ;


Para.xb = 0      ;    % x-positon of Center of Buoyancy
Para.yb = 0      ;    % y-positon of Center of Buoyancy
Para.zb = -0.02  ;    % z-positon of Center of Buoyancy

Para.rb = [Para.xb Para.yb Para.zb]' ;

%% Body positions


% Main Body Total S0;
Para.S0.r=[0.0616,0,0]'; % Position (m)

% Main Body Aft Sa;
Para.Sa.r=transpose([-(Para.Lab/4+.56) 0 0]); % Position (m)

% Main Body Middle Sm;
Para.Sm.r=transpose([.0675 0 0]); % Position (m)

% Main Body Fore Sf;
Para.Sf.r=transpose([(3*Para.R/8+.695) 0 0]); % Position (m)

% First Body S1;
Para.S1.r=transpose([.638 0 Para.La1/2+Para.R]); % Position (m)

% Second Body S2;
Para.S2.r=transpose([.371 0 Para.La2/2+Para.R]); % Position (m)

% Third Body S3;
Para.S3.r=transpose([-.453 0 Para.La3/2+Para.R]); % Position (m)

% Fourth Body S4;
Para.S4.r=transpose([-.555 0 Para.La4/2+Para.R]); % Position (m)

% Left Thruster Body Stl;
Para.Stl.r=transpose([-.59 -.17 0]); % Position (m)

% Right Thruster Body Str;
Para.Str.r=transpose([-.59 .17 0]); % Position (m)

%% Body Mass matrices


% Main Body Aft Ma;
Para.Sa.Mb = MIabC ; 

% Main Body Middle Mm;
Para.Sm.Mb = MImbC ; 

% Main Body Fore Mf;
Para.Sf.Mb = MIfbC ; 

% First Body S1;
Para.S1.Mb = MIa1C ; 

% Second Body S2;
Para.S2.Mb = MIa2C ; 

% Third Body S3;
Para.S3.Mb = MIa3C ;

% Fourth Body S4;
Para.S4.Mb = MIa4C ;

% Thruster Body TR;
Para.Stl.Mb = MIthC ;

% Thruster Body TL;
Para.Str.Mb = MIthC ;


%% Body added Mass matrices

% Main Body S0;
Para.S0.Ma = AMM ; 

% First Body S1;
Para.S1.Ma = AMM1BT ; 

% Second Body S2;
Para.S2.Ma = AMM2BT ; 

% Third Body S3;
Para.S3.Ma = AMM3BT ;

% Fourth Body S4;
Para.S4.Ma = AMM4BT ;

%% Generalized mass matrix

Para.Mg = MI_Total + AMM_Total ;


%% Generalized coriolis matrix

% Computed in RovModel.m

%% Friction matrices

% Main Body S0;
Para.S0.Kq = K_Mb ;    %Quadratic friction matrix

% First Body S1;
Para.S1.Kq = K_a1 ;    %Quadratic friction matrix

% Second Body S2;
Para.S2.Kq = K_a2 ;    %Quadratic friction matrix

% Second Body S3;
Para.S3.Kq = K_a3 ;    %Quadratic friction matrix

% Second Body S4;
Para.S4.Kq = K_a4 ;    %Quadratic friction matrix

% Left Thruster Stl;
Para.Stl.Kq = K_thl ;    %Quadratic friction matrix

% Right Thruster Str;
Para.Str.Kq = K_thr ;    %Quadratic friction matrix


%% Thruster modelling

%Thruster positions in body-fixed frame

Para.d1x = 0        ; 
Para.d1y = 0        ;
Para.d1z = 0.08     ;
Para.d2x = -0.59    ; 
Para.d2y = 0.17     ;
Para.d2z = 0        ;
Para.d3x = -0.59    ;
Para.d3y = -0.17    ;
Para.d3z = 0        ;


Para.rt1 = [Para.d1x, Para.d1y, Para.d1z]' ;
Para.rt2 = [Para.d2x, Para.d2y, Para.d2z]' ;
Para.rt3 = [Para.d3x, Para.d3y, Para.d3z]' ;


Para.rt = [Para.rt1 Para.rt2 Para.rt3] ;

%Thruster gains

Para.kt1 = 28.5  ;
Para.kt2 = 30    ;
Para.kt3 = 30    ;


%Thruster gain vectors

Para.Kt=[Para.kt1;Para.kt2;Para.kt3];

%Thruster time constants

Para.Tau1 = 0.4 ;
Para.Tau2 = 0.8 ;
Para.Tau3 = 0.8 ;


%Thruster time constant vectors

Para.Tau = [Para.Tau1;Para.Tau2;Para.Tau3] ;

% Mapping of thruster

% Para.Eb_F = [0 1 1;
%              0 0 0;
%              1 0 0];
%     
% Para.Eb_M = [0 0 0;
%              0 0 0;
%              0 -0.17 0.17]  ;
Para.Eb_F = zeros(3);

Para.Eb_M = zeros(3);

Para.Eb = [ Para.Eb_F ; Para.Eb_M ] ;

% Inverse Mapping of thruster
Para.Ebinv = pinv(Para.Eb); 

end





 
           

