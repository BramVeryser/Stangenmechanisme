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
phi = [0,0,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11]; % eerst 2 nullen voor mooie indexering
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
dphi = [0,0,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11];

ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddphi8 = zeros(size(t));
ddphi9 = zeros(size(t));
ddphi10 = zeros(size(t));
ddphi11 = zeros(size(t));
ddphi = [0,0,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10,ddphi11];
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
    PLp(k)=x(10);
    
    
    % *** velocity analysis ***
         %lus1
    A = [-BD*sin(phi3), +CD*sin(phi4),zeros(1,8);
         +BD*cos(phi3), -CD*cos(phi4),zeros(1,8);
         %lus2
         0,Ep*sin(phi4-pi/2)+CEp*sin(phi4),+EF*sin(phi5),-FpG*sin(phi6)-Fp*sin(phi6-pi/2), zeros(1,6);
         0,-Ep*cos(phi4-pi/2)-CEp*cos(phi4),-EF*cos(phi5),+FpG*cos(phi6)+Fp*cos(phi6-pi/2), zeros(1,6);
         %lus3
         zeros(1,5),-Lp8*sin(phi8-pi/2)+PLp8*sin(phi8),+Lp10*sin(phi10-pi/2)-Lp10O*sin(phi10),-OP*sin(phi11)*dphi11,-cos(phi8);
         zeros(1,5),Lp8*cos(phi8-pi/2)-PLp8*cos(phi8),-Lp10*cos(phi10-pi/2)+Lp10O*cos(phi10),+OP*cos(phi11)*dphi11,-sin(phi8);
         %lus4
         zeros(1,4),-IJ*sin(phi7)*dphi7,-Lp8*sin(phi8-pi/2)+Ip*sin(phi8-pi/2)-IpLp8*sin(phi8),+JN*sin(phi9),+Lp10*sin(phi10-pi/2)+Lp10N*sin(phi10),0,0;
         zeros(1,4),+IJ*cos(phi7)*dphi7,Lp8*cos(phi8-pi/2)-Ip*cos(phi8-pi/2)+IpLp8*cos(phi8),-JN*cos(phi9),-Lp10*cos(phi10-pi/2)-Lp10N*cos(phi10),0,0;
         %lus5
         0,+Ep*sin(phi4-pi/2)-EpK*sin(phi4),EF*sin(phi5),FpH*sin(phi6)-Fp*sin(phi6-pi/2),HI*sin(phi7),-Ip*sin(phi8-pi/2)-IpK*sin(phi8),zeros(1,4);
         0,-Ep*cos(phi4-pi/2)+EpK*cos(phi4),-EF*cos(phi5),-FpH*cos(phi6)+Fp*cos(phi6-pi/2),-HI*cos(phi7),Ip*cos(phi8-pi/2)+IpK*cos(phi8),zeros(1,4)];
     
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
    dPLp(k)=x(10);
    
    
    % *** acceleration analysis ***
    
    A = [-r3*sin(phi3(k)),  r4*sin(phi4(k));
         r3*cos(phi3(k)), -r4*cos(phi4(k))];
    B = [r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r4*cos(phi4(k))*dphi4(k)^2;
         r2*sin(phi2(k))*dphi2(k)^2-r2*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r4*sin(phi4(k))*dphi4(k)^2];
    
    x = A\B;
    % save results
    ddphi3(k) = x(1);
    ddphi4(k) = x(2);
    
    
    % *** calculate initial values for next iteration step ***
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    
    
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



