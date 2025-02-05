function f = tank60EqDT_inflow30(x,u,h)
% tank60EqDT_inflow30 - system dynamics for the uncontrolled, discrete-time 
%   version of the tank benchmark with 60 states and 30 inputs (see [1])
%
% Syntax:
%    f = tank60EqDT_inflow30(x,u,h)
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

% Authors:       Laura Luetzow, Matthias Althoff
% Written:       29-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% tanks 31-40
f(31,1) = x(31) + h*(k*sqrt(2*g)*(sqrt(x(30))-sqrt(x(31)) + u(16))); 
f(32,1) = x(32) + h*(k*sqrt(2*g)*(sqrt(x(31))-sqrt(x(32)) + u(17))); 
f(33,1) = x(33) + h*(k*sqrt(2*g)*(sqrt(x(32))-sqrt(x(33))));        
f(34,1) = x(34) + h*(k*sqrt(2*g)*(sqrt(x(33))-sqrt(x(34))));         
f(35,1) = x(35) + h*(sqrt(2*g)*(sqrt(x(34))-sqrt(x(35)))) + u(18);   
f(36,1) = x(36) + h*(k*sqrt(2*g)*(sqrt(x(35))-sqrt(x(36))));         
f(37,1) = x(37) + h*(k*sqrt(2*g)*(sqrt(x(36))-sqrt(x(37)) + u(19))); 
f(38,1) = x(38) + h*(k*sqrt(2*g)*(sqrt(x(37))-sqrt(x(38)) + u(20))); 
f(39,1) = x(39) + h*(k*sqrt(2*g)*(sqrt(x(38))-sqrt(x(39))));        
f(40,1) = x(40) + h*(k*sqrt(2*g)*(sqrt(x(39))-sqrt(x(40))));     

% tanks 41-50
f(41,1) = x(41) + h*(k*sqrt(2*g)*(sqrt(x(40))-sqrt(x(41)) + u(21))); 
f(42,1) = x(42) + h*(k*sqrt(2*g)*(sqrt(x(41))-sqrt(x(42)) + u(22))); 
f(43,1) = x(43) + h*(k*sqrt(2*g)*(sqrt(x(42))-sqrt(x(43))));        
f(44,1) = x(44) + h*(k*sqrt(2*g)*(sqrt(x(43))-sqrt(x(44))));         
f(45,1) = x(45) + h*(sqrt(2*g)*(sqrt(x(44))-sqrt(x(45)))) + u(23);   
f(46,1) = x(46) + h*(k*sqrt(2*g)*(sqrt(x(45))-sqrt(x(46))));         
f(47,1) = x(47) + h*(k*sqrt(2*g)*(sqrt(x(46))-sqrt(x(47)) + u(24))); 
f(48,1) = x(48) + h*(k*sqrt(2*g)*(sqrt(x(47))-sqrt(x(48)) + u(25))); 
f(49,1) = x(49) + h*(k*sqrt(2*g)*(sqrt(x(48))-sqrt(x(49))));        
f(50,1) = x(50) + h*(k*sqrt(2*g)*(sqrt(x(49))-sqrt(x(50))));       

% tanks 51-60
f(51,1) = x(51) + h*(k*sqrt(2*g)*(sqrt(x(50))-sqrt(x(51)) + u(26))); 
f(52,1) = x(52) + h*(k*sqrt(2*g)*(sqrt(x(51))-sqrt(x(52)) + u(27))); 
f(53,1) = x(53) + h*(k*sqrt(2*g)*(sqrt(x(52))-sqrt(x(53))));        
f(54,1) = x(54) + h*(k*sqrt(2*g)*(sqrt(x(53))-sqrt(x(54))));         
f(55,1) = x(55) + h*(sqrt(2*g)*(sqrt(x(54))-sqrt(x(55)))) + u(28);   
f(56,1) = x(56) + h*(k*sqrt(2*g)*(sqrt(x(55))-sqrt(x(56))));         
f(57,1) = x(57) + h*(k*sqrt(2*g)*(sqrt(x(56))-sqrt(x(57)) + u(29))); 
f(58,1) = x(58) + h*(k*sqrt(2*g)*(sqrt(x(57))-sqrt(x(58)) + u(30))); 
f(59,1) = x(59) + h*(k*sqrt(2*g)*(sqrt(x(58))-sqrt(x(59))));        
f(60,1) = x(60) + h*(k*sqrt(2*g)*(sqrt(x(59))-sqrt(x(60))));       

% ------------------------------ END OF CODE ------------------------------
