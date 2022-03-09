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

function [phi,dphi,ddphi] = kinematics_4bar(STANGEN,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar)
%initialisatie
AB= STANGEN(1);     BD= STANGEN(2);     CK= STANGEN(3);     Ep= STANGEN(4);
CD= STANGEN(5);     CEp= STANGEN(6);     EF= STANGEN(7);     GH= STANGEN(8);
Fp= STANGEN(9);    FpG= STANGEN(10);    HI= STANGEN(11);    IJ= STANGEN(12);
KM= STANGEN(13);    Lp8= STANGEN(14);   Ip= STANGEN(15);   KLp8= STANGEN(16);
IpK= STANGEN(17);    JN= STANGEN(18);  NO=STANGEN(19);     Lp10=STANGEN(20);
Lp10O=STANGEN(21);    OP=STANGEN(22);    ACx=STANGEN(23);    ACy=STANGEN(24);
AGx=STANGEN(25);    AGy=STANGEN(26);

Lp10N = NO - Lp10O;
IpLp8 = KLp8 - IpK;
FpH = GH - FpG;
EpK = CK - CEp;
% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
% phi1 is hoek met de wereld (keuze assenstelsel)
% phi2 is aandrijving (gekend) en heeft snelheid en versnelling
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
phi8 = zeros(size(t));
phi9 = zeros(size(t));
phi10 = zeros(size(t));
phi11 = zeros(size(t));
phi = [zeros(size(t)),zeros(size(t)),phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11]; % eerst 2 nullen voor mooie indexering
%onbekende lengte 

dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dphi8 = zeros(size(t));
dphi9 = zeros(size(t));
dphi10 = zeros(size(t));
dphi11 = zeros(size(t));
dphi = [zeros(size(t)),zeros(size(t)),dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11];

ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddphi8 = zeros(size(t));
ddphi9 = zeros(size(t));
ddphi10 = zeros(size(t));
ddphi11 = zeros(size(t));
ddphi = [zeros(size(t)),zeros(size(t)),ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10,ddphi11];
% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs',phi_init,optim_options,phi2(k),STANGEN);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    % ! volgorde nagaan !!!
    phi3(k)=x(1);
    phi4(k)=x(2);
    phi5(k)=x(3);
    phi6(k)=x(4);
    phi7(k)=x(5);
    phi8(k)=x(6);
    phi9(k)=x(7);
    phi10(k)=x(8);
    phi11(k)=x(9);
    PLp8(k)=x(10);
    
    
    % *** velocity analysis ***
         %lus1
    A = [-BD*sin(phi3(k)), +CD*sin(phi4(k)),zeros(1,8);
         +BD*cos(phi3(k)), -CD*cos(phi4(k)),zeros(1,8);
         %lus2
         0,Ep*sin(phi4(k)-pi/2)+CEp*sin(phi4(k)),+EF*sin(phi5(k)),-FpG*sin(phi6(k))-Fp*sin(phi6(k)-pi/2), zeros(1,6);
         0,-Ep*cos(phi4(k)-pi/2)-CEp*cos(phi4(k)),-EF*cos(phi5(k)),+FpG*cos(phi6(k))+Fp*cos(phi6(k)-pi/2), zeros(1,6);
         %lus3
         zeros(1,6),-Lp8*sin(phi8(k)-pi/2)+PLp8*sin(phi8(k)),+Lp10*sin(phi10(k)-pi/2)-Lp10O*sin(phi10(k)),-OP*sin(phi11(k))*dphi11(k),-cos(phi8(k));
         zeros(1,6),Lp8*cos(phi8(k)-pi/2)-PLp8*cos(phi8(k)),-Lp10*cos(phi10(k)-pi/2)+Lp10O*cos(phi10(k)),+OP*cos(phi11(k))*dphi11(k),-sin(phi8(k));
         %lus4
         zeros(1,4),-IJ*sin(phi7(k))*dphi7(k),-Lp8*sin(phi8(k)-pi/2)+Ip*sin(phi8(k)-pi/2)-IpLp8*sin(phi8(k)),+JN*sin(phi9(k)),+Lp10*sin(phi10(k)-pi/2)+Lp10N*sin(phi10(k)),0,0;
         zeros(1,4),+IJ*cos(phi7(k))*dphi7(k),Lp8*cos(phi8(k)-pi/2)-Ip*cos(phi8(k)-pi/2)+IpLp8*cos(phi8(k)),-JN*cos(phi9(k)),-Lp10*cos(phi10(k)-pi/2)-Lp10N*cos(phi10(k)),0,0;
         %lus5
         0,+Ep*sin(phi4(k)-pi/2)-EpK*sin(phi4(k)),EF*sin(phi5(k)),FpH*sin(phi6(k))-Fp*sin(phi6(k)-pi/2),HI*sin(phi7(k)),-Ip*sin(phi8(k)-pi/2)-IpK*sin(phi8(k)),zeros(1,4);
         0,-Ep*cos(phi4(k)-pi/2)+EpK*cos(phi4(k)),-EF*cos(phi5(k)),-FpH*cos(phi6(k))+Fp*cos(phi6(k)-pi/2),-HI*cos(phi7(k)),Ip*cos(phi8(k)-pi/2)+IpK*cos(phi8(k)),zeros(1,4)];
     
    B = [ AB*sin(phi2)*dphi2;
        -AB*cos(phi2)*dphi2;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0];
     
    x = A\B;
    
    % save results
    dphi3(k)=x(1);
    dphi4(k)=x(2);
    dphi5(k)=x(3);
    dphi6(k)=x(4);
    dphi7(k)=x(5);
    dphi8(k)=x(6);
    dphi9(k)=x(7);
    dphi10(k)=x(8);
    dphi11(k)=x(9);
    dPLp8(k)=x(10);
    
    
    % *** acceleration analysis ***
    %matrix A voor dynamica is A kin want d(c*dphi) => c*ddphi (matrix A) + dc*dphi (matrix B), componenten
    %afleiden naar phi moet in B staan (constanten)
    A = [-BD*sin(phi3(k)), +CD*sin(phi4(k)),zeros(1,8);
         +BD*cos(phi3(k)), -CD*cos(phi4(k)),zeros(1,8);
         %lus2
         0,Ep*sin(phi4(k)-pi/2)+CEp*sin(phi4(k)),+EF*sin(phi5(k)),-FpG*sin(phi6(k))-Fp*sin(phi6(k)-pi/2), zeros(1,6);
         0,-Ep*cos(phi4(k)-pi/2)-CEp*cos(phi4(k)),-EF*cos(phi5(k)),+FpG*cos(phi6(k))+Fp*cos(phi6(k)-pi/2), zeros(1,6);
         %lus3
         zeros(1,6),-Lp8*sin(phi8(k)-pi/2)+PLp8*sin(phi8(k)),+Lp10*sin(phi10(k)-pi/2)-Lp10O*sin(phi10(k)),-OP*sin(phi11(k))*dphi11(k),-cos(phi8(k));
         zeros(1,6),Lp8*cos(phi8(k)-pi/2)-PLp8*cos(phi8(k)),-Lp10*cos(phi10(k)-pi/2)+Lp10O*cos(phi10(k)),+OP*cos(phi11(k))*dphi11(k),-sin(phi8(k));
         %lus4
         zeros(1,4),-IJ*sin(phi7(k))*dphi7(k),-Lp8*sin(phi8(k)-pi/2)+Ip*sin(phi8(k)-pi/2)-IpLp8*sin(phi8(k)),+JN*sin(phi9(k)),+Lp10*sin(phi10(k)-pi/2)+Lp10N*sin(phi10(k)),0,0;
         zeros(1,4),+IJ*cos(phi7(k))*dphi7(k),Lp8*cos(phi8(k)-pi/2)-Ip*cos(phi8(k)-pi/2)+IpLp8*cos(phi8(k)),-JN*cos(phi9(k)),-Lp10*cos(phi10(k)-pi/2)-Lp10N*cos(phi10(k)),0,0;
         %lus5
         0,+Ep*sin(phi4(k)-pi/2)-EpK*sin(phi4(k)),EF*sin(phi5(k)),FpH*sin(phi6(k))-Fp*sin(phi6(k)-pi/2),HI*sin(phi7(k)),-Ip*sin(phi8(k)-pi/2)-IpK*sin(phi8(k)),zeros(1,4);
         0,-Ep*cos(phi4(k)-pi/2)+EpK*cos(phi4(k)),-EF*cos(phi5(k)),-FpH*cos(phi6(k))+Fp*cos(phi6(k)-pi/2),-HI*cos(phi7(k)),Ip*cos(phi8(k)-pi/2)+IpK*cos(phi8(k)),zeros(1,4)];

     
    B = [ AB*cos(phi2)*dphi2^2+AB*sin(phi2)*ddphi2-(-BD*cos(phi3(k))*dphi3(k) +CD*cos(phi4(k))*dphi4(k));
        AB*sin(phi2)*dphi2^2--AB*cos(phi2)*ddphi2-( -BD*sin(phi3(k))*dphi3(k) +CD*sin(phi4(k))*dphi4(k));
        %
        -(Ep*cos(phi4(k)-pi/2)*dphi4(k)+CEp*cos(phi4(k))*dphi4(k)+EF*cos(phi5(k))*dphi5(k)-FpG*cos(phi6(k))*dphi6(k)-Fp*cos(phi6(k)-pi/2)*dphi6(k));
        -(Ep*sin(phi4(k)-pi/2)*dphi4(k)+CEp*sin(phi4(k))*dphi4(k)+EF*sin(phi5(k))*dphi5(k)-FpG*sin(phi6(k))*dphi6(k)-Fp*sin(phi6(k)-pi/2)*dphi6(k));
        %
        -(-Lp8*cos(phi8(k)-pi/2)*dphi8(k)+PLp8*cos(phi8(k))*dphi8(k)+Lp10*cos(phi10(k)-pi/2)*dphi10(k)-Lp10O*cos(phi10(k))*dphi10(k)-OP*cos(phi11(k))*dphi11(k) +sin(phi8(k))*dphi8(k));
        -(-Lp8*sin(phi8(k)-pi/2)*dphi8(k)+PLp8*sin(phi8(k))*dphi8(k)+Lp10*sin(phi10(k)-pi/2)*dphi10(k)-Lp10O*sin(phi10(k))*dphi10(k)-OP*sin(phi11(k))*dphi11(k) -cos(phi8(k))*dphi8(k));
        %
        -(-IJ*cos(phi7(k))*dphi7(k)-Lp8*cos(phi8(k)-pi/2)*dphi8(k)+Ip*cos(phi8(k)-pi/2)*dphi8(k)-IpLp8*cos(phi8(k))*dphi8(k)+JN*cos(phi9(k))*dphi9(k)+Lp10*cos(phi10(k)-pi/2)*dphi10(k)+Lp10N*cos(phi10(k))*dphi10(k));
        -(-IJ*sin(phi7(k))*dphi7(k)-Lp8*sin(phi8(k)-pi/2)*dphi8(k)+Ip*sin(phi8(k)-pi/2)*dphi8(k)-IpLp8*sin(phi8(k))*dphi8(k)+JN*sin(phi9(k))*dphi9(k)+Lp10*sin(phi10(k)-pi/2)*dphi10(k)+Lp10N*sin(phi10(k))*dphi10(k));
        %
        -(+Ep*cos(phi4(k)-pi/2)*dphi4(k)-EpK*cos(phi4(k))*dphi4(k)+EF*cos(phi5(k))*dphi5(k)+FpH*cos(phi6(k))*dphi6(k)-Fp*cos(phi6(k)-pi/2)*dphi6(k)+HI*cos(phi7(k))*dphi7(k)-Ip*cos(phi8(k)-pi/2)*dphi8(k)-IpK*cos(phi8(k))*dphi8(k));
        -(+Ep*sin(phi4(k)-pi/2)*dphi4(k)-EpK*sin(phi4(k))*dphi4(k)+EF*sin(phi5(k))*dphi5(k)+FpH*sin(phi6(k))*dphi6(k)-Fp*sin(phi6(k)-pi/2)*dphi6(k)+HI*sin(phi7(k))*dphi7(k)-Ip*sin(phi8(k)-pi/2)*dphi8(k)-IpK*sin(phi8(k))*dphi8(k))];
    
    
    x = A\B;
    % save results
    ddphi3(k)=x(1);
    ddphi4(k)=x(2);
    ddphi5(k)=x(3);
    ddphi6(k)=x(4);
    ddphi7(k)=x(5);
    ddphi8(k)=x(6);
    ddphi9(k)=x(7);
    ddphi10(k)=x(8);
    ddphi11(k)=x(9);
    ddPLp8(k)=x(10);
    
    
    % *** calculate initial values for next iteration step ***
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
    phi8_init = phi8(k)+Ts*dphi8(k);
    phi9_init = phi9(k)+Ts*dphi9(k);
    phi10_init = phi10(k)+Ts*dphi10(k);
    phi11_init = phi11(k)+Ts*dphi11(k);
    PLp_init = PLp8(k)+Ts*dPLp8(k);
    
end % loop over positions



% *** create movie ***

% point P = fixed
P = 0;
% point S = fixed
S = r1*exp(j*phi1);
% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -1.5*r2;
y_bottom = -1.5*max(r2,r4);
x_right = r1+1.5*r4;
y_top = 1.5*max(r2,r4);

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    Q = P + r2 * exp(j*phi2(index));
    R1 = Q + r3 * exp(j*phi3(index));
    R2 = S + r4 * exp(j*phi4(index));
    
    loop1 = [P Q R1 R2 S];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_4bar
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    P = 0;
    S = r1*exp(j*phi1);
    Q = P + r2 * exp(j*phi2(index));
    R = Q + r3 * exp(j*phi3(index));
    
    figure
    assembly=[P, Q, R, S];
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(311)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(312)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(313)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(312)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(313)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(312)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
end



