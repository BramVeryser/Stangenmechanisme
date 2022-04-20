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

function [phi,dphi,ddphi] = kinematics_4bar(LINKS,phi2,dphi2,ddphi2,phi_init,t,fig_kin_4bar,Ts,S,check_kin)

%initialisation
AB = LINKS(1);        BD = LINKS(2);    CK = LINKS(3);    Ep = LINKS(4);
CD = LINKS(5);        CEp = LINKS(6);   EF = LINKS(7);    GH = LINKS(8);
Fp = LINKS(9);        FpG = LINKS(10);  HI = LINKS(11);   IJ = LINKS(12);
KM = LINKS(13);       Lp8 = LINKS(14);  Ip = LINKS(15);   KLp8 = LINKS(16);
IpK = LINKS(17);      JN = LINKS(18);   NO = LINKS(19);   Lp10 = LINKS(20);
Lp10O = LINKS(21);    OP = LINKS(22);   ACx = LINKS(23);  ACy = LINKS(24);
AGx = LINKS(25);      AGy = LINKS(26);

%extra lengts that where perviously undefined
KM = KLp8 + phi_init(10);
Lp10N = NO - Lp10O;
IpLp8 = KLp8 - IpK;
DK = CK - CD;
FpH = GH - FpG;
EpK = CK - CEp;

% allocation of the result vectors 
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
phi8 = zeros(size(t));
phi9 = zeros(size(t));
phi10 = zeros(size(t));
phi11 = zeros(size(t));
PLp8 = zeros(size(t));

dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dphi8 = zeros(size(t));
dphi9 = zeros(size(t));
dphi10 = zeros(size(t));
dphi11 = zeros(size(t));
dPLp8= zeros(size(t));

ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddphi8 = zeros(size(t));
ddphi9 = zeros(size(t));
ddphi10 = zeros(size(t));
ddphi11 = zeros(size(t));
ddPLp8= zeros(size(t));


% fsolve options
optim_options = optimset('Display','off','TolFun',1e-16);

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    [x, fval, exitflag]=fsolve('loop_closure_eqs',phi_init,optim_options,phi2(k),LINKS);
    if (exitflag ~= 1 & exitflag ~= 2 & exitflag ~= 3)
        exitflag
        display 'The fsolve exit flag was not 1, probably no convergence!'
    end
    
    % save results of fsolve
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
         %loop1
    A = [-BD*sin(phi3(k)), +CD*sin(phi4(k)),zeros(1,8);
         +BD*cos(phi3(k)), -CD*cos(phi4(k)),zeros(1,8);
         %loop2
         0,Ep*sin(phi4(k)-pi/2)+CEp*sin(phi4(k)),+EF*sin(phi5(k)),-FpG*sin(phi6(k))-Fp*sin(phi6(k)-pi/2), zeros(1,6);
         0,-Ep*cos(phi4(k)-pi/2)-CEp*cos(phi4(k)),-EF*cos(phi5(k)),+FpG*cos(phi6(k))+Fp*cos(phi6(k)-pi/2), zeros(1,6);
         %loop3
         zeros(1,5),-Lp8*sin(phi8(k)-pi/2)+PLp8(k)*sin(phi8(k)),0,+Lp10*sin(phi10(k)-pi/2)-Lp10O*sin(phi10(k)),-OP*sin(phi11(k)),-cos(phi8(k));
         zeros(1,5),Lp8*cos(phi8(k)-pi/2)-PLp8(k)*cos(phi8(k)),0,-Lp10*cos(phi10(k)-pi/2)+Lp10O*cos(phi10(k)),+OP*cos(phi11(k)),-sin(phi8(k));
         %loop4
         zeros(1,4),-IJ*sin(phi7(k)),-Lp8*sin(phi8(k)-pi/2)+Ip*sin(phi8(k)-pi/2)-IpLp8*sin(phi8(k)),+JN*sin(phi9(k)),+Lp10*sin(phi10(k)-pi/2)+Lp10N*sin(phi10(k)),0,0;
         zeros(1,4),+IJ*cos(phi7(k)),Lp8*cos(phi8(k)-pi/2)-Ip*cos(phi8(k)-pi/2)+IpLp8*cos(phi8(k)),-JN*cos(phi9(k)),-Lp10*cos(phi10(k)-pi/2)-Lp10N*cos(phi10(k)),0,0;
         %loop5
         0,+Ep*sin(phi4(k)-pi/2)-EpK*sin(phi4(k)),EF*sin(phi5(k)),FpH*sin(phi6(k))-Fp*sin(phi6(k)-pi/2),HI*sin(phi7(k)),-Ip*sin(phi8(k)-pi/2)-IpK*sin(phi8(k)),zeros(1,4);
         0,-Ep*cos(phi4(k)-pi/2)+EpK*cos(phi4(k)),-EF*cos(phi5(k)),-FpH*cos(phi6(k))+Fp*cos(phi6(k)-pi/2),-HI*cos(phi7(k)),Ip*cos(phi8(k)-pi/2)+IpK*cos(phi8(k)),zeros(1,4)];
     
    B = [ AB*sin(phi2(k))*dphi2(k);
        -AB*cos(phi2(k))*dphi2(k);
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
    
    % *** accleration analysis (A remains the same) ***
    
    B = [ AB*cos(phi2(k))*(dphi2(k))^2+AB*sin(phi2(k))*ddphi2(k)-(-BD*cos(phi3(k))*dphi3(k)^2 +CD*cos(phi4(k))*dphi4(k)^2);
          AB*sin(phi2(k))*(dphi2(k))^2-AB*cos(phi2(k))*ddphi2(k)-(-BD*sin(phi3(k))*dphi3(k)^2 +CD*sin(phi4(k))*dphi4(k)^2);
        %
        -(Ep*cos(phi4(k)-pi/2)*dphi4(k)^2+CEp*cos(phi4(k))*dphi4(k)^2+EF*cos(phi5(k))*dphi5(k)^2-FpG*cos(phi6(k))*dphi6(k)^2-Fp*cos(phi6(k)-pi/2)*dphi6(k)^2);
        -(Ep*sin(phi4(k)-pi/2)*dphi4(k)^2+CEp*sin(phi4(k))*dphi4(k)^2+EF*sin(phi5(k))*dphi5(k)^2-FpG*sin(phi6(k))*dphi6(k)^2-Fp*sin(phi6(k)-pi/2)*dphi6(k)^2);
        %
        -(-Lp8*cos(phi8(k)-pi/2)*dphi8(k)^2+PLp8(k)*cos(phi8(k))*dphi8(k)^2+Lp10*cos(phi10(k)-pi/2)*dphi10(k)^2-Lp10O*cos(phi10(k))*dphi10(k)^2-OP*cos(phi11(k))*dphi11(k)^2 +2*sin(phi8(k))*dphi8(k)*dPLp8(k));
        -(-Lp8*sin(phi8(k)-pi/2)*dphi8(k)^2+PLp8(k)*sin(phi8(k))*dphi8(k)^2+Lp10*sin(phi10(k)-pi/2)*dphi10(k)^2-Lp10O*sin(phi10(k))*dphi10(k)^2-OP*sin(phi11(k))*dphi11(k)^2 -2*cos(phi8(k))*dphi8(k)*dPLp8(k));
        %
        -(-IJ*cos(phi7(k))*dphi7(k)^2-Lp8*cos(phi8(k)-pi/2)*dphi8(k)^2+Ip*cos(phi8(k)-pi/2)*dphi8(k)^2-IpLp8*cos(phi8(k))*dphi8(k)^2+JN*cos(phi9(k))*dphi9(k)^2+Lp10*cos(phi10(k)-pi/2)*dphi10(k)^2+Lp10N*cos(phi10(k))*dphi10(k)^2);
        -(-IJ*sin(phi7(k))*dphi7(k)^2-Lp8*sin(phi8(k)-pi/2)*dphi8(k)^2+Ip*sin(phi8(k)-pi/2)*dphi8(k)^2-IpLp8*sin(phi8(k))*dphi8(k)^2+JN*sin(phi9(k))*dphi9(k)^2+Lp10*sin(phi10(k)-pi/2)*dphi10(k)^2+Lp10N*sin(phi10(k))*dphi10(k)^2);
        %
        -(+Ep*cos(phi4(k)-pi/2)*dphi4(k)^2-EpK*cos(phi4(k))*dphi4(k)^2+EF*cos(phi5(k))*dphi5(k)^2+FpH*cos(phi6(k))*dphi6(k)^2-Fp*cos(phi6(k)-pi/2)*dphi6(k)^2+HI*cos(phi7(k))*dphi7(k)^2-Ip*cos(phi8(k)-pi/2)*dphi8(k)^2-IpK*cos(phi8(k))*dphi8(k)^2);
        -(+Ep*sin(phi4(k)-pi/2)*dphi4(k)^2-EpK*sin(phi4(k))*dphi4(k)^2+EF*sin(phi5(k))*dphi5(k)^2+FpH*sin(phi6(k))*dphi6(k)^2-Fp*sin(phi6(k)-pi/2)*dphi6(k)^2+HI*sin(phi7(k))*dphi7(k)^2-Ip*sin(phi8(k)-pi/2)*dphi8(k)^2-IpK*sin(phi8(k))*dphi8(k)^2)];
    
%     rcond(A)
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
    PLp8_init = PLp8(k)+Ts*dPLp8(k);
    
    phi_init=[phi3_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,phi9_init,phi10_init,phi11_init,PLp8_init]';  
end

%phi12 is the same as phi8
phi12 = phi8;
dphi12 = dphi8;
ddphi12 = ddphi8;

%Save results in a matrix and add phi2 for later
phi = [phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11,phi12,PLp8];
dphi = [dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,dphi11,dphi12,dPLp8];
ddphi = [ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10,ddphi11,ddphi12,ddPLp8];

% *** numerical control of solutions ***

if check_kin
    position_check(phi,LINKS)
    numerical_check(Ts,t,phi,dphi,ddphi)
end

% *** create movie ***
if fig_kin_4bar 

%Ground points
A = 0;
C = ACx + j*ACy;
G = AGx + j*AGy;

% define which positions we want as frames in our movie
frames = 60;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -40*S;
y_bottom = -20*S;
x_right = 5*S;
y_top = 20*S;

figure(99)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes


% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    %Define all the points
    B = A + AB * exp(j*phi2(index));
    D = B + BD * exp(j*phi3(index));
    
    Fv = G + FpG * exp(j*phi6(index));
    F = Fv + Fp * exp(j*(phi6(index)-pi/2));
    E = F - EF * exp(j*phi5(index));
    Ev = E - Ep * exp(j*(phi4(index)-pi/2));
    
    K = Ev + (CK-CEp) * exp(j*phi4(index));
    Iv = K + IpK * exp(j*phi8(index));
    I = Iv + Ip * exp(j*(phi8(index)-pi/2));
    H = Fv + (GH-FpG) * exp(j*phi6(index));
    Lv8 = K + KLp8 * exp(j*phi8(index));
    L = Lv8 + Lp8 * exp(j*(phi8(index)-pi/2));
    Lv10 = L - Lp10 * exp(j*(phi10(index)-pi/2));
    
    J = H + (HI - IJ) * exp(j*phi7(index));
    N = J + JN * exp(j*phi9(index));
    O = N + NO * exp(j*phi10(index));
    P = O + OP * exp(j*phi11(index));
    
    M = K + KM * exp(j*phi8(index));
    
    %Define the loops
    loop1 = [A B D C];
    loop2 = [G Fv F E Ev C];
    loop3 = [Lv8 L Lv10 O P];
    loop4 = [Lv8 L Lv10 N J I Iv];
    loop5 = [Iv I H Fv F E Ev K];
    loop6 = [K M];
    
    %Plot the loops
    figure(99)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    plot(real(loop2),imag(loop2),'-o')
    plot(real(loop3),imag(loop3),'-o')
    plot(real(loop4),imag(loop4),'-o')
    plot(real(loop5),imag(loop5),'-o')
    plot(real(loop6),imag(loop6),'-o')
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film


end

% save movie
save fourbar_movie Movie
close(99)

%% *** plot figures ***
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
   
    %Define all the points
    B = A + AB * exp(j*phi2(index));
    D = B + BD * exp(j*phi3(index));
    
    Fv = G + FpG * exp(j*phi6(index));
    F = Fv + Fp * exp(j*(phi6(index)-pi/2));
    E = F - EF * exp(j*phi5(index));
    Ev = E - Ep * exp(j*(phi4(index)-pi/2));
    
    K = Ev + (CK-CEp) * exp(j*phi4(index));
    Iv = K + IpK * exp(j*phi8(index));
    I = Iv + Ip * exp(j*(phi8(index)-pi/2));
    H = Fv + (GH-FpG) * exp(j*phi6(index));
    Lv8 = K + KLp8 * exp(j*phi8(index));
    L = Lv8 + Lp8 * exp(j*(phi8(index)-pi/2));
    Lv10 = L - Lp10 * exp(j*(phi10(index)-pi/2));
    
    J = H + (HI - IJ) * exp(j*phi7(index));
    N = J + JN * exp(j*phi9(index));
    O = N + NO * exp(j*phi10(index));
    P = O + OP * exp(j*phi11(index));
    
    M = K + KM * exp(j*phi8(index));
 
    %Define the loops
    loop1 = [A B D C];
    loop2 = [G Fv F E Ev C];
    loop3 = [Lv8 L Lv10 O P];
    loop4 = [Lv8 L Lv10 N J I Iv];
    loop5 = [Iv I H Fv F E Ev K];
    loop6 = [K M];
    
    %Plot the loops
    figure()
    plot(real(loop1),imag(loop1),'ro-')
    hold on
    plot(real(loop2),imag(loop2),'ro-')
    plot(real(loop3),imag(loop3),'ro-')
    plot(real(loop4),imag(loop4),'ro-')
    plot(real(loop5),imag(loop5),'ro-')
    plot(real(loop6),imag(loop6),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    
    %plot the kinematic results
    figure
    subplot(331)
    plot(t,phi2)
    xlabel('t [s]')
    ylabel('\phi_2 [rad]')
    subplot(332)
    plot(t,phi3)
    xlabel('t [s]')
    ylabel('\phi_3 [rad]')
    subplot(333)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    xlabel('t [s]')
    
    subplot(334)
    plot(t,dphi2)
    xlabel('t [s]')
    ylabel('d\phi_2 [rad/s]')
    subplot(335)
    plot(t,dphi3)
    xlabel('t [s]')
    ylabel('d\phi_3 [rad/s]')
    subplot(336)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    
    subplot(337)
    plot(t,ddphi2)
    xlabel('t [s]')
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(338)
    plot(t,ddphi3)
    xlabel('t [s]')
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(339)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
    
    
    figure
    subplot(331)
    plot(t,phi5)
    xlabel('t [s]')
    ylabel('\phi_5 [rad]')
    subplot(332)
    plot(t,phi6)
    xlabel('t [s]')
    ylabel('\phi_6 [rad]')
    subplot(333)
    plot(t,phi7)
    ylabel('\phi_7 [rad]')
    xlabel('t [s]')
    
    
    subplot(334)
    plot(t,dphi5)
    xlabel('t [s]')
    ylabel('d\phi_5 [rad/s]')
    subplot(335)
    plot(t,dphi6)
    xlabel('t [s]')
    ylabel('d\phi_6 [rad/s]')
    subplot(336)
    plot(t,dphi7)
    ylabel('d\phi_7 [rad/s]')
    xlabel('t [s]')
    
    subplot(337)
    plot(t,ddphi5)
    xlabel('t [s]')
    ylabel('dd\phi_5 [rad/s^2]')
    subplot(338)
    plot(t,ddphi6)
    xlabel('t [s]')
    ylabel('dd\phi_6 [rad/s^2]')
    subplot(339)
    plot(t,ddphi7)
    ylabel('dd\phi_7 [rad/s^2]')
    xlabel('t [s]')
    
    figure
    subplot(331)
    plot(t,phi8)
    xlabel('t [s]')
    ylabel('\phi_8 [rad]')
    subplot(332)
    plot(t,phi9)
    xlabel('t [s]')
    ylabel('\phi_9 [rad]')
    subplot(333)
    plot(t,phi10)
    ylabel('\phi_1_0 [rad]')
    xlabel('t [s]')
    
    
    subplot(334)
    plot(t,dphi8)
    xlabel('t [s]')
    ylabel('d\phi_8 [rad/s]')
    subplot(335)
    plot(t,dphi9)
    xlabel('t [s]')
    ylabel('d\phi_9 [rad/s]')
    subplot(336)
    plot(t,dphi10)
    ylabel('d\phi_1_0 [rad/s]')
    xlabel('t [s]')
    
    
    subplot(337)
    plot(t,ddphi8)
    ylabel('dd\phi_8 [rad/s^2]')
    xlabel('t [s]')
    subplot(338)
    plot(t,ddphi9)
    ylabel('dd\phi_9 [rad/s^2]')
    xlabel('t [s]')
    subplot(339)
    plot(t,ddphi10)
    ylabel('dd\phi_1_0 [rad/s^2]')
    xlabel('t [s]')

    figure
    subplot(321)
    plot(t,phi11)
    ylabel('\phi_1_1 [rad]')
    xlabel('t [s]')
    subplot(323)
    plot(t,dphi11)
    ylabel('d\phi_1_1 [rad/s]')
    xlabel('t [s]')
    subplot(325)
    plot(t,ddphi11)
    ylabel('dd\phi_1_1 [rad/s^2]')
    xlabel('t [s]')
    
    subplot(322)
    plot(t,PLp8)
    xlabel('t [s]')
    ylabel('PL_p_8 [m]')
    subplot(324)
    plot(t,phi9)
    xlabel('t [s]')
    ylabel('dPL_p_8 [m/s]')
    subplot(326)
    plot(t,phi10)
    ylabel('ddPL_p_8 [m/s^2]')
    xlabel('t [s]')
end



