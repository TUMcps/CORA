function f = tank6EqDT_V1(x,u)
% tank6EqDT - system dynamics for the discrete-time version of the tank 
% benchmark (see Sec. VII in [1])
%
% Syntax:  
%    f = tank6EqDT(x,u,h)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - new state vector
% 
% References:
%    [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Author:        Matthias Althoff
% Written:       25-Mar-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % parameter
    k = 0.015;
    k2 = 0.01;
    g = 9.81;
    % pre-comp
%     ksqrt2g = k*sqrt(2*g);
    h = 0.05;

    % differential equations
    f(1,1) = x(1) + h*(u(1)+0.1+k2*(4-x(6))-k*sqrt(2*g)*sqrt(x(1))); % tank 1
    f(2,1) = x(2) + h*(sqrt(2*g)*(sqrt(x(1))-sqrt(x(2))));           % tank 2
    f(3,1) = x(3) + h*(k*sqrt(2*g)*(sqrt(x(2))-sqrt(x(3))));         % tank 3
    f(4,1) = x(4) + h*(k*sqrt(2*g)*(sqrt(x(3))-sqrt(x(4))));         % tank 4
    f(5,1) = x(5) + h*(k*sqrt(2*g)*(sqrt(x(4))-sqrt(x(5))));         % tank 5
    f(6,1) = x(6) + h*(k*sqrt(2*g)*(sqrt(x(5))-sqrt(x(6))));         % tank 6

%------------- END OF CODE --------------