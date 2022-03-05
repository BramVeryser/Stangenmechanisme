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
AB= STANGEN(1);     BD= STANGEN(2);     CK= STANGEN(3);     DE= STANGEN(4);
DC= STANGEN(5);     DF= STANGEN(6);     GH= STANGEN(7);     Fp= STANGEN(8);
GFp= STANGEN(9);    HI= STANGEN(10);    KM= STANGEN(11);    Lp= STANGEN(12);
Ip= STANGEN(13);    KLp= STANGEN(14);   KIp= STANGEN(15);   JN= STANGEN(16);
NO= STANGEN(17);    OP= STANGEN(18);    ACx=STANGEN(19);    ACy=STANGEN(20);
AGx=STANGEN(21);    AGy=STANGEN(22);

IpLp= KLp-KIp;
HFp = GH-GFp;
KE=   CK-DC;
% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles phi3 and phi4
phi3=phi_init(3);
phi4=phi_init(4);
phi5=phi_init(5);
phi6=phi_init(6);
phi7=phi_init(7);
phi8=phi_init(8);
phi9=phi_init(9);
phi10=phi_init(10);
phi11=phi_init(11);


%loop closure equations:
%lus 1: 1234-1
F(1)=AB*cos(phi2)+BD*cos(phi3)-DC*cos(phi4)-ACx;
F(2)=AB*sin(phi2)+BD*sin(phi3)-DC*sin(phi4)+ACy;
%lus 2: 1654-1
F(3)=GFp*cos(phi6)+Fp*cos(phi6-pi/2)+FE*cos(phi5)+ED*cos(phi5-pi/2)-CD*cos(phi4)-ACx +AGx;
F(4)=GFp*sin(phi6)+Fp*sin(phi6-pi/2)+FE*sin(phi5)+ED*sin(phi5-pi/2)-CD*sin(phi4)-ACy +AGy;
%lus 3: 81011-8 % deze bevat extra onbekende PLp: schuiverafstand
F(5)=LO*cos(phi10)+OP*cos(phi11) - PLp*cos(phi8)+Lp*cos(phi8-pi/2);
F(6)=LO*sin(phi10)+OP*sin(phi11) - PLp*sin(phi8)+Lp*sin(phi8-pi/2);
%lus 4: 81097-8
F(7)=-LN*cos(phi10)-JN*cos(phi9) + JI*cos(phi7)+ IpLp*cos(phi8)+ (Lp-Ip)*cos(phi8-pi/2);
F(8)=-LN*sin(phi10)-JN*sin(phi9) + JI*sin(phi7)+ IpLp*sin(phi8)+ (Lp-Ip)*sin(phi8-pi/2);
%lus 5: 87654-8
F(9)= HI*cos(phi7) + HFp*cos(phi6) + Fp*cos(phi6-phi/2)+ EF*cos(phi5)- KE*cos(phi4) +KIp*cos(phi8)+Ip*cos(phi8-pi/2);
F(10)=HI*sin(phi7) + HFp*sin(phi6) + Fp*sin(phi6-phi/2)+ EF*sin(phi5)- KE*sin(phi4) +KIp*sin(phi8)+Ip*sin(phi8-pi/2);
end 