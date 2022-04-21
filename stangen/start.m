%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bewegingen: 12 bar linkage, Fowler flaps
% 
%Maarten Overmeire r0797854
%Bram Veryser r0778645
%
%2021-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% program data
fig_kin_4bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1
check_kin = 1;           % Position, velocity and acceleration check
check_COG = 0;           % draw the center of gravity and the COG vectors
index = 1;               % the time step at which the assembly and COG is plotted, must be bigger than 1

%% kinematic parameters 
% bar lengths [m]
%   See pictures in report
%   p = perpendicular: e.g. Lp8 is the distance from L to link 8 (number not needed if no confusion possible)

%Scale factor 
S=0.0612;

%Bar 1(ground) coordinates of C en G with respect to A (A as origin)
AC = [-2.818; 3.031]*S;     AG=[0.376;-4.677]*S; 

%Bar 2
AB =2.645*S;

%Bar 3
BD=3.580*S;

%Bar 4
CK=8.357*S;                 Ep=1.413*S;         CD=3.215*S;     CEp=3.464*S;

%Bar 5
EF=8.242*S;

%Bar 6
GH=16.500*S;                Fp=0.578*S;         FpG=12.116*S;

%Bar 7
HI=6.797*S;                 IJ=1.3195*S; 

%Bar 8
KM=24.525*S;                Lp8=2.184*S;        Ip=4.207*S;     KLp8=17.011*S;      IpK=10.394*S;

%Bar 9
JN=6.516*S;      

%Bar 10
NO=4.9435*S;                Lp10=0.39032*S;     Lp10O=3.869*S;  

%Bar 11
OP=6.102*S; 

%Bar 12 (not included in loop equations and kinematics)
bar12 = 11*S; 

%width of the wing (for link 8 and 12) (needed for calculation of mass)
width = 25*S; 

 BARS = [AB;BD;CK;Ep;CD;CEp;EF;GH;Fp;FpG;HI;IJ;KM;Lp8;Ip;KLp8;IpK;JN;NO;Lp10;Lp10O;OP;AC;AG;bar12];

%% dynamic parameters, defined in a local frame on each of the bars.
r=0.02;       %radius of link [m]
rho = 2755;   %desity of aluminium [kg/m^3]
g = 9.81;     %gravity constant [m/s^2]

%mass of every link
m2 = AB*pi*r^2*rho;
m3 = BD*pi*r^2*rho;
m4 = (CK + Ep)*pi*r^2*rho;
m5 = EF*pi*r^2*rho;
m6 = (GH + Fp)*pi*r^2*rho;
m7 = HI*pi*r^2*rho;
m8 = (KM*0.01*width)*rho;           % A flat plate
m9 = JN*pi*r^2*rho;
m10 = (NO + Lp10)*pi*r^2*rho;
m11 = OP*pi*r^2*rho;
m12 = bar12*width*0.0075*rho;      % A flat plate

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
J12 = m12*bar12^2/12;

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
PLp8=7.514*S;           %variable slider length: initial length

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

% calculation of the kinematics (see kinematics.m)
[phi,dphi,ddphi] = kinematics(BARS,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar,S,check_kin,index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dynamics.m)
[vel,acc,F] = dynamics(phi,dphi,ddphi,phi2,dphi2,ddphi2,BARS,J,m,t,fig_dyn_4bar,S,g,check_COG);


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

%Work by the external moment in A
P = P + dphi2.*F(:,30);

%Plot P-dEkin over time, should be very small
figure()
plot(t,P-dE_kin)
ylabel('P-dE_k_i_n [W]')
xlabel('t [s]')
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
