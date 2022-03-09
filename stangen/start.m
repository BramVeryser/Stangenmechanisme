%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 0;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
AC = [-2.818; 3.031];     AG=[0.376;-4.677]; % deze zijn coordinaten tov A, al de rest zijn lengtes

AB =2.645; % stang 2

BD=3.580; % stang 3

CK=8.357;                Ep=1.413;        CD=3.215;         CEp=3.464; % stang 4

EF=8.242; % stang 5

GH=16.500;               Fp=0.578;        FpG=12.116; % stang 6

HI=6.797;       IJ=1.3195;   % stang 7

KM=24.525;               Lp8=2.184;        Ip=4.207;  KLp8=17.011;      IpK=10.394;  % stang 8

JN=6.516;  % stang 9

NO=4.9435;      Lp10=0.39032;       Lp10O=3.869;% stang 10

OP=6.102; % stang 11
% schaalfactor S
S=1;
% TODO LO, JI, GF,EF definieren, afstanden tussen scharnieren op 1 stang
STANGEN = [AB;BD;CK;Ep;CD;CEp;EF;GH;Fp;FpG;HI;IJ;KM;Lp8;Ip;KLp8;IpK;JN;NO;Lp10;Lp10O;OP;AC;AG]*S;
phi1 = 0; %grondhoek

%% dynamic parameters, defined in a local frame on each of the bars.
% zwaartepunten nog te bepalen (uitgaan van helft van de stang lijkt mij
% niet overal geldig)
% X2 = r2/2;               % X coordinates of cog (centre of gravity)
% X3 = r3/2;
% X4 = r4/2;
% 
% Y2 = 0;                  % Y coordinates of cog
% Y3 = 0.0102362;
% Y4 = 0;
% 
% m2 = r2*1.76;
% m3 = r3*1.76;
% m4 = r4*0.54;
% 
% J2 = m2*r2^2/12;
% J3 = m3*r3^2/12;
% J4 = m4*r4^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis
phi3_init = 2.67884;     % initial condition for first step of position analysis with fsolve (phi3 and phi4)
phi4_init = 3.49897;    % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi5_init = 4.04233;     %  moeten we echt een starthoek hebben voor elke
phi6_init = 3.07687;     % stang?
phi7_init = 2.58474;
phi8_init = 3.53734;
phi9_init = 3.56888;
phi10_init = 4.59245;
phi11_init = 3.31019;
PLp8=7.514*S;
phi_init=[phi3_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,phi9_init,phi10_init,phi11_init,PLp8]';

t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.3; %zelf berekenen door phi2 te meten
A = 1.019;
B = 4.043;
phi2=A*cos(omega*t+pi)+B; %tussen welke 2 hoeken ligt onze phi2?
dphi2=-omega*A*sin(omega*t);
ddphi2=-omega^2*A*cos(omega*t);

% calculation of the kinematics (see kin_4bar.m)
[phi,dphi,ddphi] = kinematics_4bar(STANGEN,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar);
%% tot hier voor eerste oefenzitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
  m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

