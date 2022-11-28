function f = SMIBdyn(x,y,u,P)
% SMIBdyn - dynamic function for a single machine infinite bus system
%           (see [1])
%
% Syntax:  
%    f = SMIBdyn(x,y,u,P)
%
% Inputs:
%    x - state vector
%    y - vector of algebraic variables
%    u - input vector
%    P - struct containing the model parameter
%
% Outputs:
%    f - time-derivative of the system state
%
% Reference:
%    [1] Ahmed El-Guindy, Control and Stability of Power Systems
%        using Reachability Analysis, PhD Thesis, TUM 2017

% Author:       Ahmed El-Guindy, Mark Wetzlinger
% Written:      22-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% State variables -------------------
delta   = x(1);
omega   = x(2);
Eq_     = x(3);

% Inputs ----------------------------
deltarl  = u(1);
omegarl  = u(2);
Eq_rl    = u(3);

% Alg. varibles ---------------------
ed      = y(1);
eq      = y(2);
id      = y(3);
iq      = y(4);
P1      = y(5);
Q1      = y(6);
Theta1  = y(7);
V1      = y(8);

x1 = (P.Tm-(P.xq-P.x_d)*id*iq)/(2*P.H*delta); 
x2 = iq; 
x3 = id/delta;
  
% Parameter Depdandant Controller ----------------------------
alpha = [ -(25*x2 - 391/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 881/266);
           (25*x2 - 391/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 615/266);
           (25*x2 - 351/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 881/266);
          -(25*x2 - 351/40)*((500*x1)/17 - 121/68)*((500*x3)/133 - 615/266);
           (25*x2 - 391/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 881/266);
          -(25*x2 - 391/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 615/266);
          -(25*x2 - 351/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 881/266);
           (25*x2 - 351/40)*((500*x1)/17 - 53/68)*((500*x3)/133 - 615/266) ];
       
K = 1.0e+04 * [0.011422432372455   1.509509521267021  -0.011003653524966;
               0.011462598370783   1.509509521189578  -0.011003653524654;
               0.011190584090320   1.447666229394653  -0.011358761159364;
               0.011230716616862   1.447664877138313  -0.011358755674191;
               0.018119719099198   1.606155453002849  -0.011331758629217;
               0.018159885100296   1.606155452990147  -0.011331758629685;
               0.017383029200531   1.567306261121143  -0.011736512340165;
               0.017423180141780   1.567308539700949  -0.011736540258299];

Kp = alpha.'*K;
vf=Kp*[delta-deltarl; omega-omegarl ;Eq_-Eq_rl]+2.4207;

% State Variables
dx(1) = P.omegaS*omega;
dx(2) = ((P.Tm-(P.xq-P.x_d)*id*iq)-iq*Eq_-P.D*omega)/(2*P.H);
dx(3) = (vf-(P.xd-P.x_d)*id-Eq_)/P.T_d0;

f = dx.';

end

%------------- END OF CODE --------------