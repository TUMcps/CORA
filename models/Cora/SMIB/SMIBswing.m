function f = SMIBswing(x,y,u,P,I)
% SMIBswing - dynamic function for a single machine infinite bus system
%             (see [1])
%
% Syntax:  
%    f = SMIBswing(x,y,u,P)
%
% Inputs:
%    x - state vector
%    y - vector of algebraic variables
%    u - input vector
%    P,I - struct containing the model parameter
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
delta1 = x(1);
omega1 = x(2);

% Alg. varibles ---------------------
P1      = y(1,:);   % active power at bus 1
Q1      = y(2,:);   % reactive power at bus 1
Theta1  = y(3,:);   % phase angle at bus 1
V1      = y(4,:);   % voltage at bus 1

Pe1 = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-Theta1));

% State Variables
dx(1) = P.omegaS*omega1;
dx(2) = (1/P.M1)*(I.Pm1-Pe1 - (P.D1*omega1));

f = dx.';

end

%------------- END OF CODE --------------