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
AB= STANGEN(1);     BD= STANGEN(2);     CK= STANGEN(3);     DE= STANGEN(4);
DC= STANGEN(5);     DF= STANGEN(6);     GH= STANGEN(7);     Fp= STANGEN(8);
GFp= STANGEN(9);    HI= STANGEN(10);    KM= STANGEN(11);    Lp= STANGEN(12);
Ip= STANGEN(13);    KLp= STANGEN(14);   KIp= STANGEN(15);   JN= STANGEN(16);
NO= STANGEN(17);    OP= STANGEN(18);    ACx=STANGEN(19);    ACy=STANGEN(20);
AGx=STANGEN(21);    AGy=STANGEN(22);

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


% loop closure equations:
F(1)=r2*cos(phi2)+r3*cos(phi3)-r4*cos(phi4)-r1*cos(phi1);
F(2)=r2*sin(phi2)+r3*sin(phi3)-r4*sin(phi4)-r1*sin(phi1);

