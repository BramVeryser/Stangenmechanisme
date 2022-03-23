function  plot_kin(t,phi1,phi2,phi3,dphi1,dphi2,dphi3,ddphi1,ddphi2,ddphi3)
%PLOT_KIN Summary of this function goes here
%   Detailed explanation goes here
figure
    subplot(331)
    plot(t,phi1)
    ylabel('x_1 [rad]')
    subplot(332)
    plot(t,phi2)
    ylabel('x_2 [rad]')
    subplot(333)
    plot(t,phi3)
    ylabel('x_3 [rad]')
    xlabel('t [s]')
    
    
    subplot(334)
    plot(t,dphi1)
    ylabel('dx_1 [rad/s]')
    subplot(335)
    plot(t,dphi2)
    ylabel('dx_2 [rad/s]')
    subplot(336)
    plot(t,dphi3)
    ylabel('dx_3 [rad/s]')
    xlabel('t [s]')
    
    
    subplot(337)
    plot(t,ddphi1)
    ylabel('ddx_1 [rad/s^2]')
    subplot(338)
    plot(t,ddphi2)
    ylabel('ddx_2 [rad/s^2]')
    subplot(339)
    plot(t,ddphi3)
    ylabel('ddx_3 [rad/s^2]')
    xlabel('t [s]')
end

