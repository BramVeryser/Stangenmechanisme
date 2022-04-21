%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bewegingen: 12 bar linkage, Fowler flaps
%NUMERICAL CHECK
% 
%Maarten Overmeire r0797854
%Bram Veryser r0778645
%
%2021-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  numerical_check(Ts,t,phi,dphi,ddphi)

%initialisation
phi3 = phi(:,2);
phi4 = phi(:,3);
phi5 = phi(:,4);
phi6 = phi(:,5);
phi7 = phi(:,6);
phi8 = phi(:,7);
phi9 = phi(:,8);
phi10 = phi(:,9);
phi11 = phi(:,10);
PLp8= phi(:,12);

dphi3 = dphi(:,2);
dphi4 = dphi(:,3);
dphi5 = dphi(:,4);
dphi6 = dphi(:,5);
dphi7 = dphi(:,6);
dphi8 = dphi(:,7);
dphi9 = dphi(:,8);
dphi10 = dphi(:,9);
dphi11 = dphi(:,10);
dPLp8= dphi(:,12);
%
ddphi3 = ddphi(:,2);
ddphi4 = ddphi(:,3);
ddphi5 = ddphi(:,4);
ddphi6 = ddphi(:,5);
ddphi7 = ddphi(:,6);
ddphi8 = ddphi(:,7);
ddphi9 = ddphi(:,8);
ddphi10= ddphi(:,9);
ddphi11= ddphi(:,10);
ddPLp8 = ddphi(:,12);

%For every variable: central defferentiation f'(x) = (f(x+1) - f(x-1)) / (2*Ts)
%We plot the exact solution, the approximation and the error (should be verry small)  
%% Phi 3

dphi3_num  = (phi3(3:size(phi3)) - phi3(1:size(phi3)-2)) / (2*Ts);
errdphi3 = dphi3_num - dphi3(2:size(phi3)-1);

ddphi3_num  = (dphi3(3:size(phi3)) - dphi3(1:size(phi3)-2)) / (2*Ts);
errddphi3 = ddphi3_num - ddphi3(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi3(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_3 [rad/s]')
    title('d\phi_3 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi3_num)
    xlabel('t [s]')
    ylabel('d\phi_3 num. [rad/s]')
    title('d\phi_3 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi3)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi3(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_3 [rad/s^2]')
    title('dd\phi_3 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi3_num)
    xlabel('t [s]]')
    ylabel('dd\phi_3 num. [rad/s^2]')
    title('dd\phi_3 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi3)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')    
%% Phi 4

dphi4_num  = (phi4(3:size(phi3)) - phi4(1:size(phi3)-2)) / (2*Ts);
errdphi4 = dphi4_num - dphi4(2:size(phi3)-1);

ddphi4_num  = (dphi4(3:size(phi3)) - dphi4(1:size(phi3)-2)) / (2*Ts);
errddphi4 = ddphi4_num - ddphi4(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi4_num)),dphi4(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_4 [rad/s]')
    title('d\phi_4 exact')
    subplot(323)
    plot(t(1:size(dphi4_num)),dphi4_num)
    xlabel('t [s]')
    ylabel('d\phi_4 num. [rad/s]')
    title('d\phi_4 numerical')
    subplot(325)
    plot(t(1:size(dphi4_num)),errdphi4)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi4(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_4 [rad/s^2]')
    title('dd\phi_4 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi4_num)
    xlabel('t [s]]')
    ylabel('dd\phi_4 num. [rad/s^2]')
    title('dd\phi_4 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi4)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% Phi 5

dphi5_num  = (phi5(3:size(phi3)) - phi5(1:size(phi3)-2)) / (2*Ts);
errdphi5 = dphi5_num - dphi5(2:size(phi3)-1);

ddphi5_num  = (dphi5(3:size(phi3)) - dphi5(1:size(phi3)-2)) / (2*Ts);
errddphi5 = ddphi5_num - ddphi5(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi5(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_5 [rad/s]')
    title('d\phi_5 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi5_num)
    xlabel('t [s]')
    ylabel('d\phi_5 num. [rad/s]')
    title('d\phi_5 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi5)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi5(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_5 [rad/s^2]')
    title('dd\phi_5 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi5_num)
    xlabel('t [s]')
    ylabel('dd\phi_5 num. [rad/s^2]')
    title('dd\phi_5 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi5)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% Phi 6

dphi6_num  = (phi6(3:size(phi3)) - phi6(1:size(phi3)-2)) / (2*Ts);
errdphi6 = dphi6_num - dphi6(2:size(phi3)-1);

ddphi6_num  = (dphi6(3:size(phi3)) - dphi6(1:size(phi3)-2)) / (2*Ts);
errddphi6 = ddphi6_num - ddphi6(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi6(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_6 [rad/s]')
    title('d\phi_6 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi6_num)
    xlabel('t [s]')
    ylabel('d\phi_6 num. [rad/s]')
    title('d\phi_6 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi6)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi6(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_6 [rad/s^2]')
    title('dd\phi_6 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi6_num)
    xlabel('t [s]')
    ylabel('dd\phi_6 num. [rad/s^2]')
    title('dd\phi_6 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi6)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% Phi 7

dphi7_num  = (phi7(3:size(phi3)) - phi7(1:size(phi3)-2)) / (2*Ts);
errdphi7 = dphi7_num - dphi7(2:size(phi3)-1);

ddphi7_num  = (dphi7(3:size(phi3)) - dphi7(1:size(phi3)-2)) / (2*Ts);
errddphi7 = ddphi7_num - ddphi7(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi7(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_7 [rad/s]')
    title('d\phi_7 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi7_num)
    xlabel('t [s]')
    ylabel('d\phi_7 num. [rad/s]')
    title('d\phi_7 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi7)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi7(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_7 [rad/s^2]')
    title('dd\phi_7 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi7_num)
    xlabel('t [s]')
    ylabel('dd\phi_7 num. [rad/s^1]')
    title('dd\phi_7 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi7)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% Phi 8 (idem Phi12)

dphi8_num  = (phi8(3:size(phi3)) - phi8(1:size(phi3)-2)) / (2*Ts);
errdphi8 = dphi8_num - dphi8(2:size(phi3)-1);

ddphi8_num  = (dphi8(3:size(phi3)) - dphi8(1:size(phi3)-2)) / (2*Ts);
errddphi8 = ddphi8_num - ddphi8(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi8(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_8 [rad/s]')
    title('d\phi_8 exact (idem \phi_1_2)')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi8_num)
    xlabel('t [s]')
    ylabel('d\phi_8 num. [rad/s]')
    title('d\phi_8 numerical (idem \phi_1_2)')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi8)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error (idem \phi_1_2)')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi8(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_8 [rad/s^2]')
    title('dd\phi_8 exact (idem \phi_1_2)')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi8_num)
    xlabel('t [s]')
    ylabel('dd\phi_8 num. [rad/s^2]')
    title('dd\phi_8 numerical (idem \phi_1_2)')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi8)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error (idem \phi_1_2)')
%% Phi 9

dphi9_num  = (phi9(3:size(phi3)) - phi9(1:size(phi3)-2)) / (2*Ts);
errdphi9 = dphi9_num - dphi9(2:size(phi3)-1);

ddphi9_num  = (dphi9(3:size(phi3)) - dphi9(1:size(phi3)-2)) / (2*Ts);
errddphi9 = ddphi9_num - ddphi9(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi9(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_9 [rad/s]')
    title('d\phi_9 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi9_num)
    xlabel('t [s]')
    ylabel('d\phi_9 num. [rad/s]')
    title('d\phi_9 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi9)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi9(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_9 [rad/s^2]')
    title('dd\phi_9 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi9_num)
    xlabel('t [s]')
    ylabel('dd\phi_9 num. [rad/s^2]')
    title('dd\phi_9 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi9)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% Phi 10

dphi10_num  = (phi10(3:size(phi3)) - phi10(1:size(phi3)-2)) / (2*Ts);
errdphi10 = dphi10_num - dphi10(2:size(phi3)-1);

ddphi10_num  = (dphi10(3:size(phi3)) - dphi10(1:size(phi3)-2)) / (2*Ts);
errddphi10 = ddphi10_num - ddphi10(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi10(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_1_0 [rad/s]')
    title('d\phi_1_0 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi10_num)
    xlabel('t [s]')
    ylabel('d\phi_1_0 num. [rad/s]')
    title('d\phi_1_0 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi10)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi10(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_1_0 [rad/s^2]')
    title('dd\phi_1_0 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi10_num)
    xlabel('t [s]')
    ylabel('dd\phi_1_0 num. [rad/s^2]')
    title('dd\phi_1_0 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi10)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% Phi 11

dphi11_num  = (phi11(3:size(phi3)) - phi11(1:size(phi3)-2)) / (2*Ts);
errdphi11 = dphi11_num - dphi11(2:size(phi3)-1);

ddphi11_num  = (dphi11(3:size(phi3)) - dphi11(1:size(phi3)-2)) / (2*Ts);
errddphi11 = ddphi11_num - ddphi11(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dphi11(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('d\phi_1_1 [rad/s]')
    title('d\phi_1_1 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dphi11_num)
    xlabel('t [s]')
    ylabel('d\phi_1_1 num. [rad/s]')
    title('d\phi_1_1 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdphi11)
    xlabel('t [s]')
    ylabel('Error [rad/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddphi11(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('dd\phi_1_1 [rad/s^2]')
    title('dd\phi_1_1 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddphi11_num)
    xlabel('t [s]')
    ylabel('dd\phi_1_1 num. [rad/s^2]')
    title('dd\phi_1_1 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddphi11)
    xlabel('t [s]')
    ylabel('Error [rad/s^2]')
    title('Error')
%% PLp8
dPLp8_num  = (PLp8(3:size(phi3)) - PLp8(1:size(phi3)-2)) / (2*Ts);
errdPLp8 = dPLp8_num - dPLp8(2:size(phi3)-1);

ddPLp8_num  = (dPLp8(3:size(phi3)) - dPLp8(1:size(phi3)-2)) / (2*Ts);
errddPLp8 = ddPLp8_num - ddPLp8(2:size(phi3)-1);

figure()
    subplot(321)
    plot(t(1:size(dphi3_num)),dPLp8(2:size(phi3,1)-1))
    xlabel('t [s]')
    ylabel('dPL_p_8 [m/s]')
    title('dPL_p_8 exact')
    subplot(323)
    plot(t(1:size(dphi3_num)),dPLp8_num)
    xlabel('t [s]')
    ylabel('dPL_p_8 num. [m/s]')
    title('dPL_p_8 numerical')
    subplot(325)
    plot(t(1:size(dphi3_num)),errdPLp8)
    xlabel('t [s]')
    ylabel('Error [m/s]')
    title('Error')
    
    subplot(322)
    plot(t(1:size(dphi3_num)),ddPLp8(2:size(phi3)-1))
    xlabel('t [s]')
    ylabel('ddPL_p_8 [m/s^2]')
    title('ddPL_p_8 exact')
    subplot(324)
    plot(t(1:size(dphi3_num)),ddPLp8_num)
    xlabel('t [s]')
    ylabel('ddPL_p_8 num. [m/s^2]')
    title('ddPL_p_8 numerical')
    subplot(326)
    plot(t(1:size(dphi3_num)),errddPLp8)
    xlabel('t [s]')
    ylabel('Error [m/s^2]')
    title('Error')
end
