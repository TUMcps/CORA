function f = tank6EqDT_inflow4(x,u,h)
% tank6EqDT - system dynamics for the uncontrolled, discrete-time version 
% of the tank benchmark with four inputs (see [1])
%
% Syntax:  
%    f = tank6EqDT(x,u,h)
%
% Inputs:
%    x - state vector
%    u - input vector
%    h - time step size
%
% Outputs:
%    f - new state vector
% 
% References:
%    [1] M. Althoff "Guaranteed State Estimation in CORA 2021", ARCH 2021

% Author:        Matthias Althoff
% Written:       25-Mar-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% parameter
k = 0.015;
g = 9.81; 

% differential equations
f(1,1) = x(1) + h*(-k*sqrt(2*g)*sqrt(x(1)) + u(1)); % tank 1
f(2,1) = x(2) + h*(sqrt(2*g)*(sqrt(x(1))-sqrt(x(2))));           % tank 2
f(3,1) = x(3) + h*(k*sqrt(2*g)*(sqrt(x(2))-sqrt(x(3))));         % tank 3
f(4,1) = x(4) + h*(k*sqrt(2*g)*(sqrt(x(3))-sqrt(x(4)) + u(2)));         % tank 4
f(5,1) = x(5) + h*(k*sqrt(2*g)*(sqrt(x(4))-sqrt(x(5)) + u(3)));         % tank 5
f(6,1) = x(6) + h*(k*sqrt(2*g)*(sqrt(x(5))-sqrt(x(6))));         % tank 6

%------------- END OF CODE --------------