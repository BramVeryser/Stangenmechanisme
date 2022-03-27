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


function [F] = ...
dynamics_4bar(phi,dphi,ddphi,phi2,dphi2,ddphi2,STANGEN,J,m,t,fig_dyn_4bar)
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
% variabelen phi dphi en ddphi moeten nog gedefinieerd worden in kin
phi3 = phi(:,1);
phi4 = phi(:,2);
phi5 = phi(:,3);
phi6 = phi(:,4);
phi7 = phi(:,5);
phi8 = phi(:,6);
phi9 = phi(:,7);
phi10 = phi(:,8);
phi11 = phi(:,9);
PLp8= phi(:,10);
%
dphi3 = dphi(:,1);
dphi4 = dphi(:,2);
dphi5 = dphi(:,3);
dphi6 = dphi(:,4);
dphi7 = dphi(:,5);
dphi8 = dphi(:,6);
dphi9 = dphi(:,7);
dphi10 = dphi(:,8);
dphi11 = dphi(:,9);
dPLp8= dphi(:,10);
%
ddphi3 = ddphi(:,1);
ddphi4 = ddphi(:,2);
ddphi5 = ddphi(:,3);
ddphi6 = ddphi(:,4);
ddphi7 = ddphi(:,5);
ddphi8 = ddphi(:,6);
ddphi9 = ddphi(:,7);
ddphi10= ddphi(:,8);
ddphi11= ddphi(:,9);
ddPLp8 = ddphi(:,10);

J2 = J(1);
J3 = J(2);
J4 = J(3);
J5 = J(4);
J6 = J(5);
J7 = J(6);
J8 = J(7);
J9 = J(8);
J10 = J(9);
J11 = J(10);

m2 = m(1);
m3 = m(2);
m4 = m(3);
m5 = m(4);
m6 = m(5);
m7 = m(6);
m8 = m(7);
m9 = m(8);
m10 = m(9);
m11 = m(10);

% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
% hefboomsarm van cog van stang i naar scharnier A...
cog2_A= -AB/2*exp(j*phi2);                                          cog2_A_x= real(cog2_A);cog2_A_y= imag(cog2_A);
cog2_B=  AB/2*exp(j*phi2);                                          cog2_B_x= real(cog2_B);cog2_B_y= imag(cog2_B);

cog3_B= -BD/2*exp(j*phi3);                                          cog3_B_x= real(cog3_B);cog3_B_y= imag(cog3_B);
cog3_D=  BD/2*exp(j*phi3);                                          cog3_D_x= real(cog3_D);cog3_D_y= imag(cog3_D);
% Zie wikipedia
cog4_C= -1/3*(CK*exp(j*phi4)+CEp*exp(j*phi4)+Ep*exp(j*(phi4-pi/2)));  cog4_C_x= real(cog4_C);cog4_C_y= imag(cog4_C);
cog4_D=  cog4_C + CD*exp(j*phi4);                                   cog4_D_x= real(cog4_D);cog4_D_y= imag(cog4_D);
cog4_E=  cog4_C + CEp*exp(j*phi4)+Ep*exp(j*(phi4-pi/2));              cog4_E_x= real(cog4_E);cog4_E_y= imag(cog4_E);
cog4_K=  cog4_C + CK*exp(j*phi4);                                   cog4_K_x= real(cog4_K);cog4_K_y= imag(cog4_K);

cog5_E= -EF/2*exp(j*phi5);                                          cog5_E_x= real(cog5_E);cog5_E_y= imag(cog5_E);
cog5_F=  EF/2*exp(j*phi5);                                          cog5_F_x= real(cog5_F);cog5_F_y= imag(cog5_F);
%3hoek
cog6_G= -1/3*(GH*exp(j*phi6)+FpG*exp(j*phi6)+Fp*exp(j*(phi6-pi/2)));  cog6_G_x= real(cog6_G);cog6_G_y= imag(cog6_G);
cog6_H=  cog6_G + GH*exp(j*phi6);                                   cog6_H_x= real(cog6_H);cog6_H_y= imag(cog6_H);
cog6_F=  cog6_G + FpG*exp(j*phi6)+Fp*exp(j*(phi6-pi/2));              cog6_F_x= real(cog6_F);cog6_F_y= imag(cog6_F);
% gn 3hoek
cog7_H= -HI/2*exp(j*phi7);                                          cog7_H_x= real(cog7_H);cog7_H_y= imag(cog7_H);
cog7_I=  HI/2*exp(j*phi7);                                          cog7_I_x= real(cog7_I);cog7_I_y= imag(cog7_I);
cog7_J=  cog7_I - IJ*exp(j*phi7) ;                                  cog7_J_x= real(cog7_J);cog7_J_y= imag(cog7_J);
%3hoek 
cog8_K= -1/3*(KM*exp(j*phi8)+IpK*exp(j*phi8)+Ip*exp(j*(phi8-pi/2)));  cog8_K_x= real(cog8_K);cog8_K_y= imag(cog8_K);
cog8_L= cog8_K + KLp8*exp(j*phi8)+Lp8*exp(j*(phi8-pi/2));              cog8_L_x= real(cog8_L);cog8_L_y= imag(cog8_L);
cog8_I= cog8_K + IpK*exp(j*phi8)+Ip*exp(j*(phi8-pi/2));               cog8_I_x= real(cog8_I);cog8_I_y= imag(cog8_I);
cog8_P= cog8_K + KLp8*exp(j*phi8)+ PLp8.*exp(j*phi8);                 cog8_P_x= real(cog8_P);cog8_P_y= imag(cog8_P);

cog9_J= -JN/2*exp(j*phi9);                                          cog9_J_x= real(cog9_J);cog9_J_y= imag(cog9_J);
cog9_N= JN/2*exp(j*phi9);                                          cog9_N_x= real(cog9_N);cog9_N_y= imag(cog9_N);
%driehoek 
cog10_N= -1/3*(NO*exp(j*phi10)+Lp10N*exp(j*phi10)+Lp10*exp(j*(phi10-pi/2)));  cog10_N_x= real(cog10_N);cog10_N_y= imag(cog10_N);
cog10_L= cog10_N + Lp10N*exp(j*phi10)+Lp10*exp(j*(phi10-pi/2));       cog10_L_x= real(cog10_L);cog10_L_y= imag(cog10_L);
cog10_O= cog10_N + NO*exp(j*phi10);                                 cog10_O_x= real(cog10_O);cog10_O_y= imag(cog10_O);

cog11_O= -OP/2*exp(j*phi11);                                        cog11_O_x= real(cog11_O);cog11_O_y= imag(cog11_O);
cog11_P=  OP/2*exp(j*phi11);                                        cog11_P_x= real(cog11_P);cog11_P_y= imag(cog11_P);


% %Controle massacentra
% index = 1;
% A = 0;
% C = ACx + j*ACy;
% G = AGx + j*AGy;
% KM = KLp8 + 7.514;
%     B = A + AB * exp(j*phi2(index));
%     D = B + BD * exp(j*phi3(index));
%     
%     Fv = G + FpG * exp(j*phi6(index));
%     F = Fv + Fp * exp(j*(phi6(index)-pi/2));
%     E = F - EF * exp(j*phi5(index));
%     Ev = E - Ep * exp(j*(phi4(index)-pi/2));
%     
%     K = Ev + (CK-CEp) * exp(j*phi4(index));
%     Iv = K + IpK * exp(j*phi8(index));
%     I = Iv + Ip * exp(j*(phi8(index)-pi/2));
%     H = Fv + (GH-FpG) * exp(j*phi6(index));
%     Lv8 = K + KLp8 * exp(j*phi8(index));
%     L = Lv8 + Lp8 * exp(j*(phi8(index)-pi/2));
%     Lv10 = L - Lp10 * exp(j*(phi10(index)-pi/2));
%     
%     J = H + (HI - IJ) * exp(j*phi7(index));
%     N = J + JN * exp(j*phi9(index));
%     O = N + NO * exp(j*phi10(index));
%     P = O + OP * exp(j*phi11(index));
%     
%     M = K + KM * exp(j*phi8(index));
% 
% figure(2)
% cog2 = A - cog2_A(index);
% stang2 = [A B];
% stang2A= [cog2 cog2+cog2_A(index)];
% stang2B= [cog2 cog2+cog2_B(index)];
% plot(real(stang2),imag(stang2),'b')
% hold on
% plot(real(stang2A),imag(stang2A),'ro-')
% plot(real(stang2B),imag(stang2B),'ro-')
% axis equal
% 
% figure(3)
% cog3 = B - cog3_B(index);
% stang3 = [B D];
% stang3B= [cog3 cog3+cog3_B(index)];
% stang3D= [cog3 cog3+cog3_D(index)];
% plot(real(stang3),imag(stang3),'b')
% hold on
% plot(real(stang3B),imag(stang3B),'ro-')
% plot(real(stang3D),imag(stang3D),'ro-')
% axis equal
% 
% figure(4)
% cog4 = C - cog4_C(index);
% stang4 = [C D K Ev E];
% stang4C= [cog4 cog4+cog4_C(index)];
% stang4D= [cog4 cog4+cog4_D(index)];
% stang4E= [cog4 cog4+cog4_E(index)];
% stang4K= [cog4 cog4+cog4_K(index)];
% plot(real(stang4),imag(stang4),'b')
% hold on
% plot(real(stang4C),imag(stang4C),'ro-')
% plot(real(stang4D),imag(stang4D),'ro-')
% plot(real(stang4E),imag(stang4E),'ro-')
% plot(real(stang4K),imag(stang4K),'ro-')
% axis equal
% 
% figure(5)
% cog5 = E - cog5_E(index);
% stang5 = [E F];
% stang5E= [cog5 cog5+cog5_E(index)];
% stang5F= [cog5 cog5+cog5_F(index)];
% plot(real(stang5),imag(stang5),'b')
% hold on
% plot(real(stang5E),imag(stang5E),'ro-')
% plot(real(stang5F),imag(stang5F),'ro-')
% axis equal
% 
% figure(6)
% cog6 = G - cog6_G(index);
% stang6 = [G H Fv F];
% stang6G= [cog6 cog6+cog6_G(index)];
% stang6H= [cog6 cog6+cog6_H(index)];
% stang6F= [cog6 cog6+cog6_F(index)];
% plot(real(stang6),imag(stang6),'b')
% hold on
% plot(real(stang6G),imag(stang6G),'ro-')
% plot(real(stang6H),imag(stang6H),'ro-')
% plot(real(stang6F),imag(stang6F),'ro-')
% axis equal
% 
% figure(7)
% cog7 = H - cog7_H(index);
% stang7 = [H I];
% stang7H= [cog7 cog7+cog7_H(index)];
% stang7I= [cog7 cog7+cog7_I(index)];
% stang7J= [cog7 cog7+cog7_J(index)];
% plot(real(stang7),imag(stang7),'b')
% hold on
% plot(real(stang7H),imag(stang7H),'ro-')
% plot(real(stang7I),imag(stang7I),'ro-')
% plot(real(stang7J),imag(stang7J),'ro-')
% axis equal
% 
% figure(8)
% cog8 = K - cog8_K(index);
% stang8 = [K M Lv8 L Lv8 Iv I];
% stang8K= [cog8 cog8+cog8_K(index)];
% stang8I= [cog8 cog8+cog8_I(index)];
% stang8L= [cog8 cog8+cog8_L(index)];
% stang8P= [cog8 cog8+cog8_P(index)];
% plot(real(stang8),imag(stang8),'b')
% hold on
% plot(real(stang8K),imag(stang8K),'ro-')
% plot(real(stang8I),imag(stang8I),'ro-')
% plot(real(stang8L),imag(stang8L),'ro-')
% plot(real(stang8P),imag(stang8P),'ro-')
% axis equal
% 
% figure(9)
% cog9 = J - cog9_J(index);
% stang9 = [J N];
% stang9J= [cog9 cog9+cog9_J(index)];
% stang9N= [cog9 cog9+cog9_N(index)];
% plot(real(stang9),imag(stang9),'b')
% hold on
% plot(real(stang9J),imag(stang9J),'ro-')
% plot(real(stang9N),imag(stang9N),'ro-')
% axis equal
% 
% figure(10)
% cog10 = N - cog10_N(index);
% stang10 = [N O Lv10 L];
% stang10N= [cog10 cog10+cog10_N(index)];
% stang10L= [cog10 cog10+cog10_L(index)];
% stang10O= [cog10 cog10+cog10_O(index)];
% plot(real(stang10),imag(stang10),'b')
% hold on
% plot(real(stang10N),imag(stang10N),'ro-')
% plot(real(stang10L),imag(stang10L),'ro-')
% plot(real(stang10O),imag(stang10O),'ro-')
% axis equal
% 
% figure(11)
% cog11 = O - cog11_O(index);
% stang11 = [O P];
% stang11O= [cog11 cog11+cog11_O(index)];
% stang11P= [cog11 cog11+cog11_P(index)];
% plot(real(stang11),imag(stang11),'b')
% hold on
% plot(real(stang11O),imag(stang11O),'ro-')
% plot(real(stang11P),imag(stang11P),'ro-')
% axis equal



% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
omega7 = [zeros(size(phi2)) zeros(size(phi2)) dphi7];
omega8 = [zeros(size(phi2)) zeros(size(phi2)) dphi8];
omega9 = [zeros(size(phi2)) zeros(size(phi2)) dphi9];
omega10 = [zeros(size(phi2)) zeros(size(phi2)) dphi10];
omega11 = [zeros(size(phi2)) zeros(size(phi2)) dphi11];

alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];
alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
alpha7 = [zeros(size(phi2)) zeros(size(phi2)) ddphi7];
alpha8 = [zeros(size(phi2)) zeros(size(phi2)) ddphi8];
alpha9 = [zeros(size(phi2)) zeros(size(phi2)) ddphi9];
alpha10 = [zeros(size(phi2)) zeros(size(phi2)) ddphi10];
alpha11 = [zeros(size(phi2)) zeros(size(phi2)) ddphi11];

% 3D model vectors
A_cog2_vec = [-cog2_A_x -cog2_A_y zeros(size(phi2))];
C_cog4_vec = [-cog4_C_x -cog4_C_y zeros(size(phi2))];
G_cog6_vec = [-cog6_G_x -cog6_G_y zeros(size(phi2))];
E_cog5_vec = [-cog5_E_x -cog5_E_y zeros(size(phi2))];
H_cog7_vec = [-cog7_H_x -cog7_H_y zeros(size(phi2))];
K_cog8_vec = [-cog8_K_x -cog8_K_y zeros(size(phi2))];
J_cog9_vec = [-cog9_J_x -cog9_J_y zeros(size(phi2))];
N_cog10_vec= [-cog10_N_x -cog10_N_y zeros(size(phi2))];
O_cog11_vec= [-cog11_O_x -cog11_O_y zeros(size(phi2))];

A_B_vec     = A_cog2_vec + [cog2_B_x cog2_B_y zeros(size(phi2))];
B_cog3_vec  = [-cog3_B_x -cog3_B_y zeros(size(phi2))];
C_E_vec     = C_cog4_vec + [cog4_E_x cog4_E_y zeros(size(phi2))];
G_H_vec     = G_cog6_vec + [cog6_H_x cog6_H_y zeros(size(phi2))];
C_K_vec     = C_cog4_vec + [cog4_K_x cog4_K_y zeros(size(phi2))];
H_J_vec     = H_cog7_vec + [cog7_J_x cog7_J_y zeros(size(phi2))];
J_N_vec     = J_cog9_vec + [cog9_N_x cog9_N_y zeros(size(phi2))];
N_O_vec     = N_cog10_vec+ [cog10_O_x cog10_O_y zeros(size(phi2))];

% acceleration vectors
acc_B =     cross(omega2,cross(omega2,A_B_vec))+cross(alpha2,A_B_vec);%
acc_E =     cross(omega4,cross(omega4,C_E_vec))+cross(alpha4,C_E_vec);%
acc_H =     cross(omega6,cross(omega6,G_H_vec))+cross(alpha6,G_H_vec);%
acc_K =     cross(omega4,cross(omega4,C_K_vec))+cross(alpha4,C_K_vec);%
acc_J =     acc_H + cross(omega7,cross(omega7,H_J_vec))+cross(alpha7,H_J_vec);%
acc_N =     acc_J + cross(omega9,cross(omega9,J_N_vec))+cross(alpha9,J_N_vec);%
acc_O =     acc_N + cross(omega10,cross(omega10,N_O_vec))+cross(alpha10,N_O_vec);%

acc_2 =       cross(omega2,cross(omega2,A_cog2_vec))+cross(alpha2,A_cog2_vec);%
acc_3 =       acc_B + cross(omega3,cross(omega3,B_cog3_vec))+cross(alpha3,B_cog3_vec);%
acc_4 =       cross(omega4,cross(omega4,C_cog4_vec))+cross(alpha4,C_cog4_vec);%
acc_5 =       acc_E + cross(omega5,cross(omega5,E_cog5_vec))+cross(alpha5,E_cog5_vec);%
acc_6 =       cross(omega6,cross(omega6,G_cog6_vec))+cross(alpha6,G_cog6_vec);%
acc_7 =       acc_H + cross(omega7,cross(omega7,H_cog7_vec))+cross(alpha7,H_cog7_vec);%
acc_8 =       acc_K + cross(omega8,cross(omega8,K_cog8_vec))+cross(alpha8,K_cog8_vec);%
acc_9 =       acc_J + cross(omega9,cross(omega9,J_cog9_vec))+cross(alpha9,J_cog9_vec);%
acc_10 =      acc_N + cross(omega10,cross(omega10,N_cog10_vec))+cross(alpha10,N_cog10_vec);%
acc_11 =      acc_O + cross(omega11,cross(omega11,O_cog11_vec))+cross(alpha11,O_cog11_vec);%

acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);
acc_5x = acc_5(:,1);
acc_5y = acc_5(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);
acc_8x = acc_8(:,1);
acc_8y = acc_8(:,2);
acc_9x = acc_9(:,1);
acc_9y = acc_9(:,2);
acc_10x = acc_10(:,1);
acc_10y = acc_10(:,2);
acc_11x = acc_11(:,1);
acc_11y = acc_11(:,2);

% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
F_A_x = zeros(size(phi2));
F_A_y = zeros(size(phi2));
F_B_x = zeros(size(phi2));
F_B_y = zeros(size(phi2));
F_C_x = zeros(size(phi2));
F_C_y = zeros(size(phi2));
F_D_x = zeros(size(phi2));
F_D_y = zeros(size(phi2));
F_E_x = zeros(size(phi2));
F_E_y = zeros(size(phi2));
F_F_x = zeros(size(phi2));
F_F_y = zeros(size(phi2));
F_G_x = zeros(size(phi2));
F_G_y = zeros(size(phi2));
F_H_x = zeros(size(phi2));
F_H_y = zeros(size(phi2));
F_I_x = zeros(size(phi2));
F_I_y = zeros(size(phi2));
F_J_x = zeros(size(phi2));
F_J_y = zeros(size(phi2));
F_K_x = zeros(size(phi2));
F_K_y = zeros(size(phi2));
F_L_x = zeros(size(phi2));
F_L_y = zeros(size(phi2));
F_N_y = zeros(size(phi2));
F_N_x = zeros(size(phi2));
F_O_y = zeros(size(phi2));
F_O_x = zeros(size(phi2));
F_P = zeros(size(phi2));
M_A = zeros(size(phi2));

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
%   A = [ 1           0            1            0            0            0            0           0           0;
%         0           1            0            1            0            0            0           0           0;
%         0           0           -1            0           -1            0            0           0           0;
%         0           0            0           -1            0           -1            0           0           0;
%         0           0            0            0            1            0            1           0           0;
%         0           0            0            0            0            1            0           1           0;
%        -cog2_P_y(k) cog2_P_x(k) -cog2_Q_y(k)  cog2_Q_x(k)  0            0            0           0           1;
%         0           0            cog3_Q_y(k) -cog3_Q_x(k)  cog3_R_y(k) -cog3_R_x(k)  0           0           0;
%         0           0            0            0           -cog4_R_y(k)  cog4_R_x(k) -cog4_S_y(k) cog4_S_x(k) 0];
%   

  A_excel = xlsread('matrix.xlsx','B2:AC21');
  A_rechts = [zeros(12,1);cos(phi8(k)-pi/2);sin(phi8(k)-pi/2);zeros(4,1);-cos(phi8(k)-pi/2);-sin(phi8(k)-pi/2)];
  A_rechts2 = [A_rechts, zeros(20,1)];
  A_boven = [A_excel,A_rechts2];
  
  M2 = [-cog2_A_y(k)    cog2_A_x(k)     -cog2_B_y(k)    cog2_B_x(k)     zeros(1,25)     1];
  M3 = [0 0 cog3_B_y(k)     -cog3_B_x(k)    0   0   cog3_D_y(k)     -cog3_D_x(k)    zeros(1,22)];
  M4 = [zeros(1,4)  -cog4_C_y(k)       cog4_C_x(k)    -cog4_D_y(k)       cog4_D_x(k)    -cog4_E_y(k)       cog4_E_x(k)    zeros(1,10)     -cog4_K_y(k)       cog4_K_x(k)    zeros(1,8)];
  M5 = [zeros(1,8)  cog5_E_y(k)       -cog5_E_x(k)    cog5_F_y(k)       -cog5_F_x(k)    zeros(1,18)];
  M6 = [zeros(1,10) -cog6_F_y(k)       cog6_F_x(k)    -cog6_G_y(k)       cog6_G_x(k)    -cog6_H_y(k)       cog6_H_x(k) zeros(1,14)];
  M7 = [zeros(1,14) cog7_H_y(k)       -cog7_H_x(k)      cog7_I_y(k)       -cog7_I_x(k)  cog7_J_y(k)       -cog7_J_x(k)  zeros(1,10)];
  M8 = [zeros(1,16)     -cog8_I_y(k)       cog8_I_x(k)  0   0   cog8_K_x(k)    -cog8_K_y(k)     -cog8_L_y(k)       cog8_L_x(k) zeros(1,4)   -cog8_P_y(k)*cos(phi8(k)-pi/2)+cog8_P_x(k)*sin(phi8(k)-pi/2)    0];
  M9 = [zeros(1,18)     -cog9_J_y(k)       cog9_J_x(k)  zeros(1,4)  -cog9_N_y(k)       cog9_N_x(k)  zeros(1,4)];
  M10 = [zeros(1,22)    cog10_L_y(k)       -cog10_L_x(k)    cog10_N_y(k)       -cog10_N_x(k)    cog10_O_y(k)       -cog10_O_x(k)    0   0];
  M11 = [zeros(1,26)    -cog11_O_y(k)   cog11_O_x(k)    cog11_P_y(k)*cos(phi8(k)-pi/2)-cog11_P_x(k)*sin(phi8(k)-pi/2) 0];
  A_onder = [M2;M3;M4;M5;M6;M7;M8;M9;M10;M11];
  
  A = [A_boven;A_onder];
  
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        m4*acc_4x(k);
        m4*acc_4y(k);
        m5*acc_5x(k);
        m5*acc_5y(k);
        m6*acc_6x(k);
        m6*acc_6y(k);
        m7*acc_7x(k);
        m7*acc_7y(k);
        m8*acc_8x(k);
        m8*acc_8y(k);
        m9*acc_9x(k);
        m9*acc_9y(k);
        m10*acc_10x(k);
        m10*acc_10y(k);
        m11*acc_11x(k);
        m11*acc_11y(k);
        J2*ddphi2(k);
        J3*ddphi3(k);
        J4*ddphi4(k)
        J5*ddphi5(k);
        J6*ddphi6(k);
        J7*ddphi7(k)
        J8*ddphi8(k);
        J9*ddphi9(k);
        J10*ddphi10(k)
        J11*ddphi11(k)];
    
    x = A\B;
    
    % save results
    F_A_x(k) = x(1);
    F_A_y(k) = x(2);
    F_B_x(k) = x(3);
    F_B_y(k) = x(4);
    F_C_x(k) = x(5);
    F_C_y(k) = x(6);
    F_D_x(k) = x(7);
    F_D_y(k) = x(8);
    F_E_x(k) = x(9);
    F_E_y(k) = x(10);
    F_F_x(k) = x(11);
    F_F_y(k) = x(12);
    F_G_x(k) = x(13);
    F_G_y(k) = x(14);
    F_H_x(k) = x(15);
    F_H_y(k) = x(16);
    F_I_x(k) = x(17);
    F_I_y(k) = x(18);
    F_J_x(k) = x(19);
    F_J_y(k) = x(20);
    F_K_x(k) = x(21);
    F_K_y(k) = x(22);
    F_L_x(k) = x(23);
    F_L_y(k) = x(24);
    F_N_x(k) = x(25);
    F_N_y(k) = x(26);
    F_O_x(k) = x(27);
    F_O_y(k) = x(28);
    F_P(k) = x(29);
    M_A(k)   = x(30);
    

end

F = [F_A_x, F_A_y, F_B_x, F_B_y, F_C_x, F_C_y, F_D_x, F_D_y, F_E_x, F_E_y, F_F_x, F_F_y, F_G_x, F_G_y, F_H_x, F_H_y, F_I_x, F_I_y, F_J_x,  F_J_y, F_K_x, F_K_y, F_L_x, F_L_y, F_N_y, F_N_x, F_O_y, F_O_x, F_P, M_A];

% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar
    
    figure
    subplot(221)
    plot(F_A_x,F_A_y),grid
    xlabel('F_A_x [N]')
    ylabel('F_A_y [N]')
    axis tight
    subplot(222)
    plot(F_B_x,F_B_y),grid
    xlabel('F_B_x [N]')
    ylabel('F_B_y [N]')
    axis tight
    subplot(223)
    plot(F_C_x,F_C_y),grid
    xlabel('F_C_x [N]')
    ylabel('F_C_y [N]')
    axis tight
    subplot(224)
    plot(F_D_x,F_D_y),grid
    xlabel('F_D_x [N]')
    ylabel('F_D_y [N]')
    axis tight
    
    figure
    plot(t,M_A)
    ylabel('M_A [N-m]')
    xlabel('t [s]')
    
end


