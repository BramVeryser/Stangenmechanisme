%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bewegingen: 12 bar linkage, Fowler flaps
%DYNAMICS
% 
%Maarten Overmeire r0797854
%Bram Veryser r0778645
%
%2021-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vel,acc,F] = ...
dynamics_4bar(phi,dphi,ddphi,phi2,dphi2,ddphi2,BARS,J,m,t,fig_dyn_4bar,S,g,check_COG)


%initialisation
AB= BARS(1);     BD= BARS(2);     CK= BARS(3);     Ep= BARS(4);
CD= BARS(5);     CEp= BARS(6);     EF= BARS(7);     GH= BARS(8);
Fp= BARS(9);    FpG= BARS(10);    HI= BARS(11);    IJ= BARS(12);
KM= BARS(13);    Lp8= BARS(14);   Ip= BARS(15);   KLp8= BARS(16);
IpK= BARS(17);    JN= BARS(18);  NO=BARS(19);     Lp10=BARS(20);
Lp10O=BARS(21);    OP=BARS(22);    ACx=BARS(23);    ACy=BARS(24);
AGx=BARS(25);    AGy=BARS(26);    bar12 = BARS(27);

%extra lengts that where perviously undefined
Lp10N = NO - Lp10O;
IpLp8 = KLp8 - IpK;
FpH = GH - FpG;
EpK = CK - CEp;

%more initialisation
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
J12 = J(11);

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
m12 = m(11);

%% COG
% cogi_P = vector from the centre of gravity of bar i to point P
% The method used is based on the same principle as is used in kinematics for the movie

cog2_A= -AB/2*exp(j*phi2);                                          cog2_A_x= real(cog2_A);cog2_A_y= imag(cog2_A);
cog2_B=  AB/2*exp(j*phi2);                                          cog2_B_x= real(cog2_B);cog2_B_y= imag(cog2_B);

cog3_B= -BD/2*exp(j*phi3);                                          cog3_B_x= real(cog3_B);cog3_B_y= imag(cog3_B);
cog3_D=  BD/2*exp(j*phi3);                                          cog3_D_x= real(cog3_D);cog3_D_y= imag(cog3_D);

% Bar 4 is assumed to be a triangle (see report)
cog4_C= -1/3*(CK*exp(j*phi4)+CEp*exp(j*phi4)+Ep*exp(j*(phi4-pi/2)));  cog4_C_x= real(cog4_C);cog4_C_y= imag(cog4_C);
cog4_D=  cog4_C + CD*exp(j*phi4);                                   cog4_D_x= real(cog4_D);cog4_D_y= imag(cog4_D);
cog4_E=  cog4_C + CEp*exp(j*phi4)+Ep*exp(j*(phi4-pi/2));              cog4_E_x= real(cog4_E);cog4_E_y= imag(cog4_E);
cog4_K=  cog4_C + CK*exp(j*phi4);                                   cog4_K_x= real(cog4_K);cog4_K_y= imag(cog4_K);

cog5_E= -EF/2*exp(j*phi5);                                          cog5_E_x= real(cog5_E);cog5_E_y= imag(cog5_E);
cog5_F=  EF/2*exp(j*phi5);                                          cog5_F_x= real(cog5_F);cog5_F_y= imag(cog5_F);

% Bar 6 is assumed to be a triangle (see report)
cog6_G= -1/3*(GH*exp(j*phi6)+FpG*exp(j*phi6)+Fp*exp(j*(phi6-pi/2)));  cog6_G_x= real(cog6_G);cog6_G_y= imag(cog6_G);
cog6_H=  cog6_G + GH*exp(j*phi6);                                   cog6_H_x= real(cog6_H);cog6_H_y= imag(cog6_H);
cog6_F=  cog6_G + FpG*exp(j*phi6)+Fp*exp(j*(phi6-pi/2));              cog6_F_x= real(cog6_F);cog6_F_y= imag(cog6_F);

cog7_H= -HI/2*exp(j*phi7);                                          cog7_H_x= real(cog7_H);cog7_H_y= imag(cog7_H);
cog7_I=  HI/2*exp(j*phi7);                                          cog7_I_x= real(cog7_I);cog7_I_y= imag(cog7_I);
cog7_J=  cog7_I - IJ*exp(j*phi7) ;                                  cog7_J_x= real(cog7_J);cog7_J_y= imag(cog7_J);

% Bar 8 is assumed to be a triangle (see report)
cog8_K= -1/3*(KM*exp(j*phi8)+IpK*exp(j*phi8)+Ip*exp(j*(phi8-pi/2)));  cog8_K_x= real(cog8_K);cog8_K_y= imag(cog8_K);
cog8_L= cog8_K + KLp8*exp(j*phi8)+Lp8*exp(j*(phi8-pi/2));              cog8_L_x= real(cog8_L);cog8_L_y= imag(cog8_L);
cog8_I= cog8_K + IpK*exp(j*phi8)+Ip*exp(j*(phi8-pi/2));               cog8_I_x= real(cog8_I);cog8_I_y= imag(cog8_I);
cog8_P= cog8_K + KLp8*exp(j*phi8)+ PLp8.*exp(j*phi8);                 cog8_P_x= real(cog8_P);cog8_P_y= imag(cog8_P);
cog8_M = cog8_K + KM*exp(j*phi8);                                      cog8_M_x= real(cog8_M);cog8_M_y= imag(cog8_M);      

cog9_J= -JN/2*exp(j*phi9);                                          cog9_J_x= real(cog9_J);cog9_J_y= imag(cog9_J);
cog9_N= JN/2*exp(j*phi9);                                          cog9_N_x= real(cog9_N);cog9_N_y= imag(cog9_N);

% Bar 10 is assumed to be a triangle (see report)
cog10_N= -1/3*(NO*exp(j*phi10)+Lp10N*exp(j*phi10)+Lp10*exp(j*(phi10-pi/2)));  cog10_N_x= real(cog10_N);cog10_N_y= imag(cog10_N);
cog10_L= cog10_N + Lp10N*exp(j*phi10)+Lp10*exp(j*(phi10-pi/2));       cog10_L_x= real(cog10_L);cog10_L_y= imag(cog10_L);
cog10_O= cog10_N + NO*exp(j*phi10);                                 cog10_O_x= real(cog10_O);cog10_O_y= imag(cog10_O);

cog11_O= -OP/2*exp(j*phi11);                                        cog11_O_x= real(cog11_O);cog11_O_y= imag(cog11_O);
cog11_P=  OP/2*exp(j*phi11);                                        cog11_P_x= real(cog11_P);cog11_P_y= imag(cog11_P);

cog12_P= -bar12/2*exp(j*phi8);                                    cog12_P_x= real(cog12_P);cog12_P_y= imag(cog12_P);

%% Check center of gravity
%visualy check the COG by plotting the bars and the COGi_. vectors
if check_COG
    index = 1;
    A = 0;
    C = ACx + j*ACy;
    G = AGx + j*AGy;
    KM = KLp8 + 7.514*S;
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
        Q = P + bar12 * exp(j*phi8(index));
        
    %Bar 2
    figure()
    cog2 = A - cog2_A(index);
    stang2 = [A B];
    stang2A= [cog2 cog2+cog2_A(index)];
    stang2B= [cog2 cog2+cog2_B(index)];
    plot(real(stang2),imag(stang2),'b')
    hold on
    plot(real(stang2A),imag(stang2A),'ro-')
    plot(real(stang2B),imag(stang2B),'ro-')
    axis equal
    
    %Bar 3
    figure()
    cog3 = B - cog3_B(index);
    stang3 = [B D];
    stang3B= [cog3 cog3+cog3_B(index)];
    stang3D= [cog3 cog3+cog3_D(index)];
    plot(real(stang3),imag(stang3),'b')
    hold on
    plot(real(stang3B),imag(stang3B),'ro-')
    plot(real(stang3D),imag(stang3D),'ro-')
    axis equal
    
    %Bar 4
    figure()
    cog4 = C - cog4_C(index);
    stang4 = [C D K Ev E];
    stang4C= [cog4 cog4+cog4_C(index)];
    stang4D= [cog4 cog4+cog4_D(index)];
    stang4E= [cog4 cog4+cog4_E(index)];
    stang4K= [cog4 cog4+cog4_K(index)];
    plot(real(stang4),imag(stang4),'b')
    hold on
    plot(real(stang4C),imag(stang4C),'ro-')
    plot(real(stang4D),imag(stang4D),'ro-')
    plot(real(stang4E),imag(stang4E),'ro-')
    plot(real(stang4K),imag(stang4K),'ro-')
    axis equal
    
    %Bar 5
    figure()
    cog5 = E - cog5_E(index);
    stang5 = [E F];
    stang5E= [cog5 cog5+cog5_E(index)];
    stang5F= [cog5 cog5+cog5_F(index)];
    plot(real(stang5),imag(stang5),'b')
    hold on
    plot(real(stang5E),imag(stang5E),'ro-')
    plot(real(stang5F),imag(stang5F),'ro-')
    axis equal
   
    %Bar 6
    figure()
    cog6 = G - cog6_G(index);
    stang6 = [G H Fv F];
    stang6G= [cog6 cog6+cog6_G(index)];
    stang6H= [cog6 cog6+cog6_H(index)];
    stang6F= [cog6 cog6+cog6_F(index)];
    plot(real(stang6),imag(stang6),'b')
    hold on
    plot(real(stang6G),imag(stang6G),'ro-')
    plot(real(stang6H),imag(stang6H),'ro-')
    plot(real(stang6F),imag(stang6F),'ro-')
    axis equal
    
    %Bar 7
    figure()
    cog7 = H - cog7_H(index);
    stang7 = [H I];
    stang7H= [cog7 cog7+cog7_H(index)];
    stang7I= [cog7 cog7+cog7_I(index)];
    stang7J= [cog7 cog7+cog7_J(index)];
    plot(real(stang7),imag(stang7),'b')
    hold on
    plot(real(stang7H),imag(stang7H),'ro-')
    plot(real(stang7I),imag(stang7I),'ro-')
    plot(real(stang7J),imag(stang7J),'ro-')
    axis equal

    %Bar 8
    figure()
    cog8 = K - cog8_K(index);
    stang8 = [K M Lv8 L Lv8 Iv I];
    stang8K= [cog8 cog8+cog8_K(index)];
    stang8I= [cog8 cog8+cog8_I(index)];
    stang8L= [cog8 cog8+cog8_L(index)];
    stang8P= [cog8 cog8+cog8_P(index)];
    plot(real(stang8),imag(stang8),'b')
    hold on
    plot(real(stang8K),imag(stang8K),'ro-')
    plot(real(stang8I),imag(stang8I),'ro-')
    plot(real(stang8L),imag(stang8L),'ro-')
    plot(real(stang8P),imag(stang8P),'ro-')
    axis equal
    
    %Bar 9
    figure()
    cog9 = J - cog9_J(index);
    stang9 = [J N];
    stang9J= [cog9 cog9+cog9_J(index)];
    stang9N= [cog9 cog9+cog9_N(index)];
    plot(real(stang9),imag(stang9),'b')
    hold on
    plot(real(stang9J),imag(stang9J),'ro-')
    plot(real(stang9N),imag(stang9N),'ro-')
    axis equal
    
    %Bar 10
    figure()
    cog10 = N - cog10_N(index);
    stang10 = [N O Lv10 L];
    stang10N= [cog10 cog10+cog10_N(index)];
    stang10L= [cog10 cog10+cog10_L(index)];
    stang10O= [cog10 cog10+cog10_O(index)];
    plot(real(stang10),imag(stang10),'b')
    hold on
    plot(real(stang10N),imag(stang10N),'ro-')
    plot(real(stang10L),imag(stang10L),'ro-')
    plot(real(stang10O),imag(stang10O),'ro-')
    axis equal
    
    %Bar 11
    figure()
    cog11 = O - cog11_O(index);
    stang11 = [O P];
    stang11O= [cog11 cog11+cog11_O(index)];
    stang11P= [cog11 cog11+cog11_P(index)];
    plot(real(stang11),imag(stang11),'b')
    hold on
    plot(real(stang11O),imag(stang11O),'ro-')
    plot(real(stang11P),imag(stang11P),'ro-')
    axis equal
 
    %Bar 12
    figure()
    cog12 = P - cog12_P(index);
    stang12 = [P Q];
    stang12P= [cog12 cog12+cog12_P(index)];
    plot(real(stang12),imag(stang12),'b')
    hold on
    plot(real(stang12P),imag(stang12P),'ro-')
    axis equal
    
end

%% 3D omega (dphi) and alpha (ddphi) vectors)
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
omega12 = omega8;

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
alpha12 = alpha8;

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
P_cog12_vec= [-cog12_P_x -cog12_P_y zeros(size(phi2))];

A_B_vec     = A_cog2_vec + [cog2_B_x cog2_B_y zeros(size(phi2))];
B_cog3_vec  = [-cog3_B_x -cog3_B_y zeros(size(phi2))];
C_E_vec     = C_cog4_vec + [cog4_E_x cog4_E_y zeros(size(phi2))];
G_H_vec     = G_cog6_vec + [cog6_H_x cog6_H_y zeros(size(phi2))];
C_K_vec     = C_cog4_vec + [cog4_K_x cog4_K_y zeros(size(phi2))];
H_J_vec     = H_cog7_vec + [cog7_J_x cog7_J_y zeros(size(phi2))];
J_N_vec     = J_cog9_vec + [cog9_N_x cog9_N_y zeros(size(phi2))];
N_O_vec     = N_cog10_vec+ [cog10_O_x cog10_O_y zeros(size(phi2))];
K_P_vec     = [(KLp8+PLp8).*cos(phi8)    (KLp8+PLp8).*sin(phi8) zeros(size(phi2))];
K_M_vec     = [(KM)*cos(phi8)    (KM)*sin(phi8) zeros(size(phi2))];
O_P_vec     = [OP.*cos(phi11)   OP.*sin(phi11)  zeros(size(phi2))];
O_P_vec     = [OP.*cos(phi11)   OP.*sin(phi11)  zeros(size(phi2))];

% velocity vectors
vel_B = cross(omega2,A_B_vec);
vel_E = cross(omega4,C_E_vec);
vel_H = cross(omega6,G_H_vec);
vel_K = cross(omega4,C_K_vec);
vel_J = vel_H + cross(omega7,H_J_vec);
vel_N = vel_J + cross(omega9,J_N_vec);
vel_O = vel_N + cross(omega10,N_O_vec);
vel_P = vel_K + cross(omega8,K_P_vec);
vel_P2 = vel_O + cross(omega11,O_P_vec);
vel_M = vel_K + cross(omega8,K_M_vec);

vel_2 =     cross(omega2,A_cog2_vec);
vel_3 =     vel_B + cross(omega3,B_cog3_vec);
vel_4 =     cross(omega4,C_cog4_vec);
vel_5 =     vel_E + + cross(omega5,E_cog5_vec);
vel_6 =     cross(omega6,G_cog6_vec);
vel_7 =     vel_H + cross(omega7,H_cog7_vec);
vel_8 =     vel_K + cross(omega8,K_cog8_vec);
vel_9 =     vel_J + cross(omega9,J_cog9_vec);
vel_10 =    vel_N + cross(omega10,N_cog10_vec);
vel_11 =    vel_O + cross(omega11,O_cog11_vec);
vel_12 =    vel_P2 + cross(omega12,P_cog12_vec);

vel = [vel_2, vel_3, vel_4,vel_5, vel_6, vel_7,vel_8, vel_9, vel_10,vel_11,vel_12];

% acceleration vectors
acc_B =     cross(omega2,cross(omega2,A_B_vec))+cross(alpha2,A_B_vec);%
acc_E =     cross(omega4,cross(omega4,C_E_vec))+cross(alpha4,C_E_vec);%
acc_H =     cross(omega6,cross(omega6,G_H_vec))+cross(alpha6,G_H_vec);%
acc_K =     cross(omega4,cross(omega4,C_K_vec))+cross(alpha4,C_K_vec);%
acc_J =     acc_H + cross(omega7,cross(omega7,H_J_vec))+cross(alpha7,H_J_vec);%
acc_N =     acc_J + cross(omega9,cross(omega9,J_N_vec))+cross(alpha9,J_N_vec);%
acc_O =     acc_N + cross(omega10,cross(omega10,N_O_vec))+cross(alpha10,N_O_vec);%
acc_P =     acc_O + cross(omega11,cross(omega11,O_P_vec))+cross(alpha11,O_P_vec);

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
acc_12 =      acc_P + cross(omega12,cross(omega12,P_cog12_vec))+cross(alpha12,P_cog12_vec);

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
acc_12x = acc_12(:,1);
acc_12y = acc_12(:,2);

acc = [acc_2,acc_3,acc_4,acc_5,acc_6,acc_7,acc_8,acc_9,acc_10,acc_11,acc_12];

%% allocate matrices for force (F) and moment (M)
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
F_P_p = zeros(size(phi2));
M_A = zeros(size(phi2));
F_P_y = zeros(size(phi2));
F_P_x = zeros(size(phi2));
M_P = zeros(size(phi2));

%% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
  A = [ 1              0             1              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              1             0              1              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             -1             0              0              0             -1             0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              -1             0              0             0              -1             0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              1              0             1              0              1              0              0              0              0              0             0              0              0              0              0              0              1              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              1             0              1              0              1              0              0              0              0             0              0              0              0              0              0              0              1              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              -1             0              -1             0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              -1             0              -1             0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              1              0              1              0             1              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              1              0              1             0              1              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             -1             0              -1             0              -1             0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              -1             0              -1             0              -1             0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              1              0              0              0              -1             0              1              0               0              0               0               0               cos(phi8(k)-pi/2)                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              1              0              0              0              -1             0              1               0              0               0               0               sin(phi8(k)-pi/2)                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              1              0              0              0              0              0               1              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              1              0              0              0              0               0              1               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              -1             0               -1             0               -1              0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              -1              0              -1              0               -1              0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               1               0               0                                                               0   1               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               1               0                                                               0   0               1               0;
        -cog2_A_y(k)   cog2_A_x(k)   -cog2_B_y(k)   cog2_B_x(k)    0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               1   0               0               0;
        0              0             cog3_B_y(k)    -cog3_B_x(k)   0              0             cog3_D_y(k)    -cog3_D_x(k)   0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0; 
        0              0             0              0              -cog4_C_y(k)   cog4_C_x(k)   -cog4_D_y(k)   cog4_D_x(k)    -cog4_E_y(k)   cog4_E_x(k)    0              0              0              0             0              0              0              0              0              0              -cog4_K_y(k)   cog4_K_x(k)    0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              cog5_E_y(k)    -cog5_E_x(k)   cog5_F_y(k)    -cog5_F_x(k)   0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              -cog6_F_y(k)   cog6_F_x(k)    -cog6_G_y(k)   cog6_G_x(k)   -cog6_H_y(k)   cog6_H_x(k)    0              0              0              0              0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             cog7_H_y(k)    -cog7_H_x(k)   cog7_I_y(k)    -cog7_I_x(k)   cog7_J_y(k)    -cog7_J_x(k)   0              0              0              0               0              0               0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              -cog8_I_y(k)   cog8_I_x(k)    0              0              cog8_K_y(k)    -cog8_K_x(k)   -cog8_L_y(k)   cog8_L_x(k)     0              0               0               0               -cog8_P_y(k)*cos(phi8(k)-pi/2)+cog8_P_x(k)*sin(phi8(k)-pi/2)    0   0               0               1;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              -cog9_J_y(k)   cog9_J_x(k)    0              0              0              0               -cog9_N_y(k)   cog9_N_x(k)     0               0               0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              cog10_L_y(k)   -cog10_L_x(k)   cog10_N_y(k)   -cog10_N_x(k)   cog10_O_y(k)    -cog10_O_x(k)   0                                                               0   0               0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               -cog11_O_y(k)   cog11_O_x(k)    0                                                               0   -cog11_P_y(k)   cog11_P_x(k)    0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               -cos(phi8(k)-pi/2)                                              0   -1              0               0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               -sin(phi8(k)-pi/2)                                              0   0               -1              0;
        0              0             0              0              0              0             0              0              0              0              0              0              0              0             0              0              0              0              0              0              0              0              0              0               0              0               0               0               cog12_P_y(k)*cos(phi8(k)-pi/2)-cog12_P_x(k)*sin(phi8(k)-pi/2)   0   cog12_P_y(k)    -cog12_P_x(k)   -1
      ];
  
  
  B = [ m2*acc_2x(k);
        m2*(acc_2y(k)+g);
        m3*acc_3x(k);
        m3*(acc_3y(k)+g);
        m4*acc_4x(k);
        m4*(acc_4y(k)+g);
        m5*acc_5x(k);
        m5*(acc_5y(k)+g);
        m6*acc_6x(k);
        m6*(acc_6y(k)+g);
        m7*acc_7x(k);
        m7*(acc_7y(k)+g);
        m8*acc_8x(k);
        m8*(acc_8y(k)+g);
        m9*acc_9x(k);
        m9*(acc_9y(k)+g);
        m10*acc_10x(k);
        m10*(acc_10y(k)+g);
        m11*acc_11x(k);
        m11*(acc_11y(k)+g);
        J2*ddphi2(k);
        J3*ddphi3(k);
        J4*ddphi4(k);
        J5*ddphi5(k);
        J6*ddphi6(k);
        J7*ddphi7(k);
        J8*ddphi8(k);
        J9*ddphi9(k);
        J10*ddphi10(k);
        J11*ddphi11(k);
        m12*acc_12x(k);
        m12*(acc_12y(k)+g);
        J12*ddphi8(k)];
    

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
    F_P_p(k) = x(29);
    M_A(k)   = x(30);
    F_P_x(k) = x(31);
    F_P_y(k) = x(32);
    M_P(k) = x(33);
    

end

F = [F_A_x, F_A_y, F_B_x, F_B_y, F_C_x, F_C_y, F_D_x, F_D_y, F_E_x, F_E_y, F_F_x, F_F_y, F_G_x, F_G_y, F_H_x, F_H_y, F_I_x, F_I_y, F_J_x,  F_J_y, F_K_x, F_K_y, F_L_x, F_L_y, F_N_y, F_N_x, F_O_y, F_O_x, F_P_p, M_A];

% **********************
%% ** plot figures ***
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
    subplot(221)
    plot(F_E_x,F_E_y),grid
    xlabel('F_E_x [N]')
    ylabel('F_E_y [N]')
    axis tight
    subplot(222)
    plot(F_F_x,F_F_y),grid
    xlabel('F_F_x [N]')
    ylabel('F_F_y [N]')
    axis tight
    subplot(223)
    plot(F_G_x,F_G_y),grid
    xlabel('F_G_x [N]')
    ylabel('F_G_y [N]')
    axis tight
    subplot(224)
    plot(F_H_x,F_H_y),grid
    xlabel('F_H_x [N]')
    ylabel('F_H_y [N]')
    axis tight
    
    figure
    subplot(221)
    plot(F_I_x,F_I_y),grid
    xlabel('F_I_x [N]')
    ylabel('F_I_y [N]')
    axis tight
    subplot(222)
    plot(F_J_x,F_J_y),grid
    xlabel('F_J_x [N]')
    ylabel('F_J_y [N]')
    axis tight
    subplot(223)
    plot(F_K_x,F_K_y),grid
    xlabel('F_K_x [N]')
    ylabel('F_K_y [N]')
    axis tight
    subplot(224)
    plot(F_L_x,F_L_y),grid
    xlabel('F_L_x [N]')
    ylabel('F_L_y [N]')
    axis tight
    
    figure
    subplot(221)
    plot(F_N_x,F_N_y),grid
    xlabel('F_N_x [N]')
    ylabel('F_N_y [N]')
    axis tight
    subplot(222)
    plot(F_O_x,F_O_y),grid
    xlabel('F_O_x [N]')
    ylabel('F_O_y [N]')
    axis tight    
    subplot(223)
    plot(F_P_p.*cos(phi8-pi/2),F_P_p.*sin(phi8-pi/2)),grid
    xlabel('F_P_p [N]')
    ylabel('F_P_p [N]')
    axis tight
    subplot(224)
    plot(F_P_x,F_P_y),grid
    xlabel('F_P_x [N]')
    ylabel('F_P_y [N]')
    axis tight
    
    figure
    plot(t,M_P)
    ylabel('M_P [Nm]')
    xlabel('t [s]')
    
    figure
    plot(t,M_A)
    ylabel('M_A [Nm]')
    xlabel('t [s]')
    
end


