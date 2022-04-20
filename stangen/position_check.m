function position_check(phi,LINKS)
%POSITION_CHECK Summary of this function goes here
%   Detailed explanation goes here
%initialisation
AB = LINKS(1);        BD = LINKS(2);    CK = LINKS(3);    Ep = LINKS(4);
CD = LINKS(5);        CEp = LINKS(6);   EF = LINKS(7);    GH = LINKS(8);
Fp = LINKS(9);        FpG = LINKS(10);  HI = LINKS(11);   IJ = LINKS(12);
KM = LINKS(13);       Lp8 = LINKS(14);  Ip = LINKS(15);   KLp8 = LINKS(16);
IpK = LINKS(17);      JN = LINKS(18);   NO = LINKS(19);   Lp10 = LINKS(20);
Lp10O = LINKS(21);    OP = LINKS(22);   ACx = LINKS(23);  ACy = LINKS(24);
AGx = LINKS(25);      AGy = LINKS(26);

%extra lengts that where perviously undefined
Lp10N = NO - Lp10O;
IpLp8 = KLp8 - IpK;
DK = CK - CD;
FpH = GH - FpG;
EpK = CK - CEp;
HJ = HI - IJ;

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

plot(real(Error))