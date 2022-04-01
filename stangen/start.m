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
fig_kin_4bar = 0;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

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

m2 = 1;
m3 = 1;
m4 = 1;
m5 = 1;
m6 = 1;
m7 = 1;
m8 = 1;
m9 = 1;
m10 = 1;
m11 = 1;
m = [m2 m3 m4 m5 m6 m7 m8 m9 m10 m11];

J2 = 0.001;
J3 = 0.001;
J4 = 0.001;
J5 = 0.001;
J6 = 0.001;
J7 = 0.001;
J8 = 0.001;
J9 = 0.001;
J10 = 0.001;
J11 = 0.001;

J = [J2 J3 J4 J5 J6 J7 J8 J9 J10 J11];





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
t_end = 15;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.2; %zelf berekenen door phi2 te meten
A = 1.019;
B = 4.043;
phi2=A*cos(omega*t+pi)+B; %tussen welke 2 hoeken ligt onze phi2?
dphi2=-omega*A*sin(omega*t+pi);
ddphi2=-omega^2*A*cos(omega*t+pi);

% calculation of the kinematics (see kin_4bar.m)
[phi,dphi,ddphi] = kinematics_4bar(STANGEN,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar,Ts);
%% tot hier voor eerste oefenzitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[F] = dynamics_4bar(phi,dphi,ddphi,phi2,dphi2,ddphi2,STANGEN,J,m,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2.5 Controle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variation of kinetic energy 
% dE_kin/dt = P
% P: vermogen geleverd door de motor P=T*omega
% del(mv^2/2) 
% Methode van virtuele arbeid 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

