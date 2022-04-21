function position_check(phi,BARS,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bewegingen: 12 bar linkage, Fowler flaps
%POSITION CHECK
% 
%Maarten Overmeire r0797854
%Bram Veryser r0778645
%
%2021-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialisation
AB = BARS(1);        BD = BARS(2);    CK = BARS(3);    Ep = BARS(4);
CD = BARS(5);        CEp = BARS(6);   EF = BARS(7);    GH = BARS(8);
Fp = BARS(9);        FpG = BARS(10);  HI = BARS(11);   IJ = BARS(12);
KM = BARS(13);       Lp8 = BARS(14);  Ip = BARS(15);   KLp8 = BARS(16);
IpK = BARS(17);      JN = BARS(18);   NO = BARS(19);   Lp10 = BARS(20);
Lp10O = BARS(21);    OP = BARS(22);   ACx = BARS(23);  ACy = BARS(24);
AGx = BARS(25);      AGy = BARS(26);

%extra lengts that where perviously undefined
Lp10N = NO - Lp10O;
IpLp8 = KLp8 - IpK;
DK = CK - CD;
FpH = GH - FpG;
EpK = CK - CEp;
HJ = HI - IJ;

%more initialisation
phi2 = phi(:,1);
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


%Two ways of calculation point P: via bar: 2 - 3 - 4 - 8 and 6- 7 - 9 - 10 - 11
%Ground points
A = 0;
C = ACx + j*ACy;
G = AGx + j*AGy;
P_1 = zeros(size(phi2));
P_2 = zeros(size(phi2));
for i = 1:size(phi2)
    
    %2 - 3 - 4 - 8
    B = A + AB * exp(j*phi2(i));
    D = B + BD * exp(j*phi3(i));
    K = D + DK * exp(j*phi4(i));
    P_1(i) = K + (KLp8+PLp8(i)) * exp(j*phi8(i));
    
    %6- 7 - 9 - 10 - 11
    H = G + GH * exp(j*phi6(i));
    J = H + HJ * exp(j*phi7(i));
    N = J + JN * exp(j*phi9(i));
    O = N + NO * exp(j*phi10(i));
    P_2(i) = O + OP * exp(j*phi11(i)); 

end

Error = P_1 - P_2;
figure()
plot(t,real(Error)) %check the x-difference over time
xlabel('t [s]')
ylabel('P_1-P_2 [m]')
title("Position check")