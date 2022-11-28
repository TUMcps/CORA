function f = tank30EqDT_inflow15(x,u,h)
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
f(1,1) = x(1) + h*(-k*sqrt(2*g)*sqrt(x(1)) + u(1));                  % tank 1
f(2,1) = x(2) + h*(sqrt(2*g)*(sqrt(x(1))-sqrt(x(2))));               % tank 2
f(3,1) = x(3) + h*(k*sqrt(2*g)*(sqrt(x(2))-sqrt(x(3))));             % tank 3
f(4,1) = x(4) + h*(k*sqrt(2*g)*(sqrt(x(3))-sqrt(x(4)) + u(2)));      % tank 4
f(5,1) = x(5) + h*(k*sqrt(2*g)*(sqrt(x(4))-sqrt(x(5)) + u(3)));      % tank 5
f(6,1) = x(6) + h*(k*sqrt(2*g)*(sqrt(x(5))-sqrt(x(6))));             % tank 6
f(7,1) = x(7) + h*(sqrt(2*g)*(sqrt(x(6))-sqrt(x(7)))) + u(4);        % tank 7
f(8,1) = x(8) + h*(k*sqrt(2*g)*(sqrt(x(7))-sqrt(x(8))));             % tank 8
f(9,1) = x(9) + h*(k*sqrt(2*g)*(sqrt(x(8))-sqrt(x(9)) + u(5)));      % tank 9
f(10,1) = x(10) + h*(k*sqrt(2*g)*(sqrt(x(9))-sqrt(x(10)) + u(6)));   % tank 10
f(11,1) = x(11) + h*(k*sqrt(2*g)*(sqrt(x(10))-sqrt(x(11))));         % tank 11
f(12,1) = x(12) + h*(k*sqrt(2*g)*(sqrt(x(11))-sqrt(x(12))));         % tank 12
f(13,1) = x(13) + h*(sqrt(2*g)*(sqrt(x(12))-sqrt(x(13)))) + u(7);    % tank 13
f(14,1) = x(14) + h*(k*sqrt(2*g)*(sqrt(x(13))-sqrt(x(14))));         % tank 14
f(15,1) = x(15) + h*(k*sqrt(2*g)*(sqrt(x(14))-sqrt(x(15)) + u(8)));  % tank 15
f(16,1) = x(16) + h*(k*sqrt(2*g)*(sqrt(x(15))-sqrt(x(16)) + u(9)));  % tank 16
f(17,1) = x(17) + h*(k*sqrt(2*g)*(sqrt(x(16))-sqrt(x(17))));         % tank 17
f(18,1) = x(18) + h*(k*sqrt(2*g)*(sqrt(x(17))-sqrt(x(18))));         % tank 18
f(19,1) = x(19) + h*(sqrt(2*g)*(sqrt(x(18))-sqrt(x(19)))) + u(10);   % tank 19
f(20,1) = x(20) + h*(k*sqrt(2*g)*(sqrt(x(19))-sqrt(x(20))));         % tank 20
f(21,1) = x(21) + h*(k*sqrt(2*g)*(sqrt(x(20))-sqrt(x(21)) + u(11))); % tank 21
f(22,1) = x(22) + h*(k*sqrt(2*g)*(sqrt(x(21))-sqrt(x(22)) + u(12))); % tank 22
f(23,1) = x(23) + h*(k*sqrt(2*g)*(sqrt(x(22))-sqrt(x(23))));         % tank 23
f(24,1) = x(24) + h*(k*sqrt(2*g)*(sqrt(x(23))-sqrt(x(24))));         % tank 24
f(25,1) = x(25) + h*(sqrt(2*g)*(sqrt(x(24))-sqrt(x(25)))) + u(13);   % tank 25
f(26,1) = x(26) + h*(k*sqrt(2*g)*(sqrt(x(25))-sqrt(x(26))));         % tank 26
f(27,1) = x(27) + h*(k*sqrt(2*g)*(sqrt(x(26))-sqrt(x(27)) + u(14))); % tank 27
f(28,1) = x(28) + h*(k*sqrt(2*g)*(sqrt(x(27))-sqrt(x(28)) + u(15))); % tank 28
f(29,1) = x(29) + h*(k*sqrt(2*g)*(sqrt(x(28))-sqrt(x(29))));         % tank 29
f(30,1) = x(30) + h*(k*sqrt(2*g)*(sqrt(x(29))-sqrt(x(30))));         % tank 30

%------------- END OF CODE --------------