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
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
S=0.0612; %Schaalfactor

AC = [-2.818; 3.031]*S;     AG=[0.376;-4.677]*S; % deze zijn coordinaten tov A, al de rest zijn lengtes

AB =2.645*S; % stang 2

BD=3.580*S; % stang 3

CK=8.357*S;                Ep=1.413*S;        CD=3.215*S;         CEp=3.464*S; % stang 4

EF=8.242*S; % stang 5

GH=16.500*S;               Fp=0.578*S;        FpG=12.116*S; % stang 6

HI=6.797*S;       IJ=1.3195*S;   % stang 7

KM=24.525*S;               Lp8=2.184*S;        Ip=4.207*S;  KLp8=17.011*S;      IpK=10.394*S;  % stang 8

JN=6.516*S;  % stang 9

NO=4.9435*S;      Lp10=0.39032*S;       Lp10O=3.869*S;% stang 10

OP=6.102*S; % stang 11

STANGEN = [AB;BD;CK;Ep;CD;CEp;EF;GH;Fp;FpG;HI;IJ;KM;Lp8;Ip;KLp8;IpK;JN;NO;Lp10;Lp10O;OP;AC;AG];
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
r=0.02;
rho = 2755;   %aluminium kg/m^3

m2 = AB*pi*r^2*rho;
m3 = BD*pi*r^2*rho;
m4 = (CK + Ep)*pi*r^2*rho;
m5 = EF*pi*r^2*rho;
m6 = (GH + Fp)*pi*r^2*rho;
m7 = HI*pi*r^2*rho;
m8 = (KM+Lp8+Ip)*pi*r^2*rho;
m9 = JN*pi*r^2*rho;
m10 = (NO + Lp10)*pi*r^2*rho;
m11 = OP*pi*r^2*rho;
m = [m2 m3 m4 m5 m6 m7 m8 m9 m10 m11];

J2 = m2*r^2/12;
J3 = m3*r^2/12;
J4 = m4*r^2/12;
J5 = m5*r^2/12;
J6 = m6*r^2/12;
J7 = m7*r^2/12;
J8 = m8*r^2/12;
J9 = m9*r^2/12;
J10 = m10*r^2/12;
J11 = m11*r^2/12;

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
t_end = 30;                    % end time of simulation
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
[phi,dphi,ddphi] = kinematics_4bar(STANGEN,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar,Ts,S);
%% tot hier voor eerste oefenzitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[vel,acc,F] = dynamics_4bar(phi,dphi,ddphi,phi2,dphi2,ddphi2,STANGEN,J,m,t,fig_dyn_4bar,S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2.5 Controle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variation of kinetic energy 
% dE_kin/dt = P
% P: vermogen geleverd door de motor P=T*omega
% del(mv^2/2) 
% Methode van virtuele arbeid
vel_norm=zeros(size(dphi));acc_norm=zeros(size(phi));
dE_kin =zeros(size(vel,1),1);dE_kin0=zeros(size(vel,1),1);
for k = 1:10 %itereer over alle stangen: x en y componenten van de vel en acc (2x10)
    vel_norm(:,k) = sqrt(vel(:,2*k-1).^2+vel(:,2*k).^2);acc_norm(:,k) = sqrt(acc(:,2*k-1).^2+acc(:,2*k).^2);
    dE_kin0 = dE_kin0 + m(ceil(k/2))*vel(:,k).*acc(:,k)+J(ceil(k/2))*dphi(:,ceil(k/2)).*ddphi(:,ceil(k/2));
    dE_kin = dE_kin + m(k)*vel_norm(:,k).*acc_norm(:,k)+J(k)*dphi(k).*ddphi(:,k); 
end
P = dphi2.*F(:,30);
figure()
plot(P-dE_kin)
hold on
plot(P-dE_kin0)
title("controle P-dE_{kin} over tijd")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

