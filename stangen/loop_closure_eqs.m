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


function F=loop_closure_eqs(phi_init,phi2,STANGEN)
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
%IpLp= KLp-KIp;
%HFp = GH-GFp;
%KE=   CK-CD;
% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles phi3 and phi4
phi3=phi_init(1);
phi4=phi_init(2);
phi5=phi_init(3);
phi6=phi_init(4);
phi7=phi_init(5);
phi8=phi_init(6);
phi9=phi_init(7);
phi10=phi_init(8);
phi11=phi_init(9);
PLp8=phi_init(10);

%loop closure equations:
%lus 1: 1234-1
F(1)=AB*cos(phi2)+BD*cos(phi3)-CD*cos(phi4)-ACx;
F(2)=AB*sin(phi2)+BD*sin(phi3)-CD*sin(phi4)-ACy;

%lus 2: 1654-1
F(3)=FpG*cos(phi6)+Fp*cos(phi6-pi/2)-EF*cos(phi5)-Ep*cos(phi4-pi/2)-CEp*cos(phi4)-ACx +AGx; %ED->Ep
F(4)=FpG*sin(phi6)+Fp*sin(phi6-pi/2)-EF*sin(phi5)-Ep*sin(phi4-pi/2)-CEp*sin(phi4)-ACy +AGy;

%lus 3: 81011-8 % deze bevat extra onbekende PLp: schuiverafstand
F(5)=Lp8*cos(phi8-pi/2)-Lp10*cos(phi10-pi/2)+Lp10O*cos(phi10)+OP*cos(phi11)-PLp8*cos(phi8);
F(6)=Lp8*sin(phi8-pi/2)-Lp10*sin(phi10-pi/2)+Lp10O*sin(phi10)+OP*sin(phi11)-PLp8*sin(phi8);
%F(5)=LO*cos(phi10)+OP*cos(phi11) - PLp*cos(phi8)+Lp*cos(phi8-pi/2);% extra Lp8 en Lp10
%F(6)=LO*sin(phi10)+OP*sin(phi11) - PLp*sin(phi8)+Lp*sin(phi8-pi/2);

%lus 4: 81097-8
F(7)=Lp8*cos(phi8-pi/2)-Lp10*cos(phi10-pi/2)-Lp10N*cos(phi10)-JN*cos(phi9)+IJ*cos(phi7)-Ip*cos(phi8-pi/2)+IpLp8*cos(phi8);
F(8)=Lp8*sin(phi8-pi/2)-Lp10*sin(phi10-pi/2)-Lp10N*sin(phi10)-JN*sin(phi9)+IJ*sin(phi7)-Ip*sin(phi8-pi/2)+IpLp8*sin(phi8);
%F(7)=-LN*cos(phi10)-JN*cos(phi9) + JI*cos(phi7)+ IpLp*cos(phi8)+ (Lp-Ip)*cos(phi8-pi/2);
%F(8)=-LN*sin(phi10)-JN*sin(phi9) + JI*sin(phi7)+ IpLp*sin(phi8)+ (Lp-Ip)*sin(phi8-pi/2);

%lus 5: 87654-8
F(9)=Ip*cos(phi8-pi/2)-HI*cos(phi7)-FpH*cos(phi6)+Fp*cos(phi6-pi/2)-EF*cos(phi5)-Ep*cos(phi4-pi/2)+EpK*cos(phi4)+IpK*cos(phi8);
F(10)=Ip*sin(phi8-pi/2)-HI*sin(phi7)-FpH*sin(phi6)+Fp*sin(phi6-pi/2)-EF*sin(phi5)-Ep*sin(phi4-pi/2)+EpK*sin(phi4)+IpK*sin(phi8);
%F(9)= -HI*cos(phi7) - HFp*cos(phi6) + Fp*cos(phi6-phi/2)- EF*cos(phi5)- KE*cos(phi4) +KIp*cos(phi8)+Ip*cos(phi8-pi/2);
%F(10)=-HI*sin(phi7) - HFp*sin(phi6) + Fp*sin(phi6-phi/2)- EF*sin(phi5)- KE*sin(phi4) +KIp*sin(phi8)+Ip*sin(phi8-pi/2);
end 