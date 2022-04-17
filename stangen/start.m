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
num_control = 0;


% kinematic parameters (link lengths)
%   See pictures 
%   p = perpendicular: e.g. Lp8 is the distance from L to link 8 (number not needed if no confusion possible)

%Scale factor 
S=0.0612;

%Link 1(ground) coordinates of C en G with A as origin
AC = [-2.818; 3.031]*S;     AG=[0.376;-4.677]*S; 

%Link 2
AB =2.645*S;

 % link 3
BD=3.580*S;

%Link 4
CK=8.357*S;                 Ep=1.413*S;         CD=3.215*S;     CEp=3.464*S;

%Link 5
EF=8.242*S;

%Link 6
GH=16.500*S;                Fp=0.578*S;         FpG=12.116*S;

%Link 7
HI=6.797*S;                 IJ=1.3195*S; 

%Link 8
KM=24.525*S;                Lp8=2.184*S;        Ip=4.207*S;     KLp8=17.011*S;      IpK=10.394*S;

%Link 9
JN=6.516*S;      

%Link 10
NO=4.9435*S;                Lp10=0.39032*S;     Lp10O=3.869*S;  

% link 11
OP=6.102*S; 

%Link 12 (not included in loop equations and kinematics)
link12 = 11*S; 

%width of the wing (for link 8 and 12) (needed of calculation of mass)
width = 25*S; 

 LINKS = [AB;BD;CK;Ep;CD;CEp;EF;GH;Fp;FpG;HI;IJ;KM;Lp8;Ip;KLp8;IpK;JN;NO;Lp10;Lp10O;OP;AC;AG;link12];

%% dynamic parameters, defined in a local frame on each of the bars.
r=0.02;       %radius of link [m]
rho = 2755;   %desity of aluminium [kg/m^3]
g = 9.81;     %gravity constant

%mass of every link
m2 = AB*pi*r^2*rho;
m3 = BD*pi*r^2*rho;
m4 = (CK + Ep)*pi*r^2*rho;
m5 = EF*pi*r^2*rho;
m6 = (GH + Fp)*pi*r^2*rho;
m7 = HI*pi*r^2*rho;
m8 = (KM*0.01*width)*rho;
m9 = JN*pi*r^2*rho;
m10 = (NO + Lp10)*pi*r^2*rho;
m11 = OP*pi*r^2*rho;
m12 = link12*width*0.0075*rho;
m = [m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12];


%moment of inertia for each link
J2 = m2*AB^2/12;
J3 = m3*BD^2/12;
J4 = m4*CK^2/12;
J5 = m5*EF^2/12;
J6 = m6*GH^2/12;
J7 = m7*HI^2/12;
J8 = m8*KM^2/12;
J9 = m9*JN^2/12;
J10 = m10*NO^2/12;
J11 = m11*OP^2/12;
J12 = m12*link12^2/12;

J = [J2 J3 J4 J5 J6 J7 J8 J9 J10 J11 J12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis, inital values for every angle 
phi3_init = 2.67884;    
phi4_init = 3.49897;     
phi5_init = 4.04233;     
phi6_init = 3.07687;     
phi7_init = 2.58474;
phi8_init = 3.53734;
phi9_init = 3.56888;
phi10_init = 4.59245;
phi11_init = 3.31019;
PLp8=7.514*S;           %variable slider length

phi_init=[phi3_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,phi9_init,phi10_init,phi11_init,PLp8]';

%Time parameters 
t_begin = 0;                   % start time of simulation
t_end = 15;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = pi/15; 
A = 1.019;
B = 4.043;
phi2=A*cos(omega*t+pi)+B; 
dphi2=-omega*A*sin(omega*t+pi);
ddphi2=-omega^2*A*cos(omega*t+pi);

% calculation of the kinematics (see kinematics_4bar.m)
[phi,dphi,ddphi] = kinematics_4bar(LINKS,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar,Ts,S,num_control);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dynamics_4bar.m)
[vel,acc,F] = dynamics_4bar(phi,dphi,ddphi,phi2,dphi2,ddphi2,LINKS,J,m,t,fig_dyn_4bar,S,g);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2.5 Controle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variation of kinetic energy 
dE_kin =zeros(size(vel,1),1);
P = zeros(size(vel,1),1);

for k = 1:11 %itteration over every link 
    vel_dot_acc = zeros(size(vel(:,1)));
    for h = 1:size(vel(:,1))                                        %dot product of velocity and acceleration for every time step
        dot_h = dot(vel(h,3*k-2:3*k),acc(h,3*k-2:3*k));
        vel_dot_acc(h) = dot_h;
    end
    dE_kin = dE_kin + m(k)*vel_dot_acc+J(k)*dphi(:,k).*ddphi(:,k);
    P = P + (-ones(size(vel,1),1)*m(k)*g).*vel(:,3*k-1);            %gravity 
end

%work by the external moment in A
P = P + dphi2.*F(:,30);

figure()
plot(P-dE_kin)
hold on
title("Check P-dE_{kin} over time")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig_kin_4bar
    figure
    load fourbar_movie Movie
    movie(Movie)
end
