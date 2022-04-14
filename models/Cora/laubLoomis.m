function dx = laubLoomis(x,u)
% laubLoomis - dynamic equation for the Laub-Loomis benchmark 
%              (see Sec. 3.2 in [1])
%
% Syntax:  
%    dx = laubLoomis(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%    [1] F. Immler, “ARCH-COMP19 Category Report: Continuous and Hybrid 
%        Systems with Nonlinear Dynamics", 2019

% Author:        Niklas Kochdumper
% Written:       19-June-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    dx(1,1) = 1.4*x(3)-0.9*x(1);
    dx(2,1) = 2.5*x(5)-1.5*x(2);
    dx(3,1) = 0.6*x(7)-0.8*x(3)*x(2);
    dx(4,1) = 2.0-1.3*x(4)*x(3);
    dx(5,1) = 0.7*x(1)-1.0*x(4)*x(5);
    dx(6,1) = 0.3*x(1)-3.1*x(6);
    dx(7,1) = 1.8*x(6)-1.5*x(7)*x(2);

%------------- END OF CODE --------------