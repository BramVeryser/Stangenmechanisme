%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bewegingen: 12 bar linkage, Fowler flaps
%LOOP CLOSURE EQUATIONS
% 
%Maarten Overmeire r0797854
%Bram Veryser r0778645
%
%2021-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=loop_closure_eqs(phi_init,phi2,BARS)
%initialisation
AB= BARS(1);     BD= BARS(2);     CK= BARS(3);     Ep= BARS(4);
CD= BARS(5);     CEp= BARS(6);     EF= BARS(7);     GH= BARS(8);
Fp= BARS(9);    FpG= BARS(10);    HI= BARS(11);    IJ= BARS(12);
KM= BARS(13);    Lp8= BARS(14);   Ip= BARS(15);   KLp8= BARS(16);
IpK= BARS(17);    JN= BARS(18);  NO=BARS(19);     Lp10=BARS(20);
Lp10O=BARS(21);    OP=BARS(22);    ACx=BARS(23);    ACy=BARS(24);
AGx=BARS(25);    AGy=BARS(26);

%extra lengts that where perviously undefined
Lp10N = NO - Lp10O;
IpLp8 = KLp8 - IpK;
FpH = GH - FpG;
EpK = CK - CEp;


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
%loop 1: 1234-1
F(1)=AB*cos(phi2)+BD*cos(phi3)-CD*cos(phi4)-ACx;
F(2)=AB*sin(phi2)+BD*sin(phi3)-CD*sin(phi4)-ACy;

%loop 2: 1654-1
F(3)=FpG*cos(phi6)+Fp*cos(phi6-pi/2)-EF*cos(phi5)-Ep*cos(phi4-pi/2)-CEp*cos(phi4)-ACx +AGx; %ED->Ep
F(4)=FpG*sin(phi6)+Fp*sin(phi6-pi/2)-EF*sin(phi5)-Ep*sin(phi4-pi/2)-CEp*sin(phi4)-ACy +AGy;

%loop 3: 81011-8 % deze bevat extra onbekende PLp: schuiverafstand
F(5)=Lp8*cos(phi8-pi/2)-Lp10*cos(phi10-pi/2)+Lp10O*cos(phi10)+OP*cos(phi11)-PLp8*cos(phi8);
F(6)=Lp8*sin(phi8-pi/2)-Lp10*sin(phi10-pi/2)+Lp10O*sin(phi10)+OP*sin(phi11)-PLp8*sin(phi8);

%loop 4: 81097-8
F(7)=Lp8*cos(phi8-pi/2)-Lp10*cos(phi10-pi/2)-Lp10N*cos(phi10)-JN*cos(phi9)+IJ*cos(phi7)-Ip*cos(phi8-pi/2)+IpLp8*cos(phi8);
F(8)=Lp8*sin(phi8-pi/2)-Lp10*sin(phi10-pi/2)-Lp10N*sin(phi10)-JN*sin(phi9)+IJ*sin(phi7)-Ip*sin(phi8-pi/2)+IpLp8*sin(phi8);

%loop 5: 87654-8
F(9)=Ip*cos(phi8-pi/2)-HI*cos(phi7)-FpH*cos(phi6)+Fp*cos(phi6-pi/2)-EF*cos(phi5)-Ep*cos(phi4-pi/2)+EpK*cos(phi4)+IpK*cos(phi8);
F(10)=Ip*sin(phi8-pi/2)-HI*sin(phi7)-FpH*sin(phi6)+Fp*sin(phi6-pi/2)-EF*sin(phi5)-Ep*sin(phi4-pi/2)+EpK*sin(phi4)+IpK*sin(phi8);
end 