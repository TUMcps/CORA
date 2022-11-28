function f = SMIBcon(x,y,u,P)
% SMIBcon - constraint function for a single machine infinite bus system
%           (see [1])
%
% Syntax:  
%    f = SMIBcon(x,y,u,P)
%
% Inputs:
%    x - state vector
%    y - vector of algebraic variables
%    u - input vector
%    P - struct containing the model parameter
%
% Outputs:
%    f - vector storing the constraint viaolation
%
% Reference:
%    [1] Ahmed El-Guindy, Control and Stability of Power Systems
%        using Reachability Analysis, PhD Thesis, TUM 2017

% Author:       Ahmed El-Guindy, Mark Wetzlinger
% Written:      22-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if strcmp(P.mode,'normal')
    xs = P.xs;
elseif strcmp(P.mode,'fault')
    xs = P.xs3;
end

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

% Alg. Variables;
dy(1) = -ed+V1*sin(delta-Theta1);
dy(2) = -eq+V1*cos(delta-Theta1);
dy(3) = -ed+P.xq*iq-P.ra*id;
dy(4) = -eq+(Eq_)-P.ra*iq-P.x_d*id;
dy(5) = -P1+ed*id+eq*iq;
dy(6) = -Q1+eq*id-ed*iq;
dy(7) = V1*P.v2/xs*sin(Theta1-P.theta2)-P1;
dy(8) = V1*V1/xs-V1*P.v2/xs*cos(Theta1-P.theta2)-Q1;

f = dy.';

end

%------------- END OF CODE --------------