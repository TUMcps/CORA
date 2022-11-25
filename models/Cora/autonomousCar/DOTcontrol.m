function [carInput,k] = DOTcontrol(x,u,y)
% DOTcontrol - provides the steering angle and acceleration of the vehicle
%              (see [1])
%
% Syntax:  
%    f = DOTcontrol(x,u,y)
%
% Inputs:
%    x - state vector
%    u - reference trajectory
%    y - sensor noise
%
% Outputs:
%    carInput - control inputs for the car
%    k - control parameter
%
% References:
%    [1] M. Althoff and J. M. Dolan. Online verification of automated
%        road vehicles using reachability analysis.
%        IEEE Transactions on Robotics, 30(4):903-918, 2014. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DOTBicycleDynamics_controlled_SRX_velEq

% Author:       Matthias Althoff
% Written:      23-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%control parameters
k(1) = 0.2;
k(2) = 2;
%k(3) = 2;
k(3) = 0.3;
k(4) = 1;
k(5) = 10; 


%compute steering angle and acceleration input from reference trajectory
% steering angle delta
delta = k(1)*(-sin(u(3))*(u(1) - (x(5) + y(1))) + cos(u(3))*(u(2) - (x(6) + y(2)))) ...
    + k(2)*(u(3) - (x(2) + y(3))) + k(3)*(u(4) - (x(3) + y(4)));
% longitudinal acceleration ax 
ax = k(4)*(cos(u(3))*(u(1) - (x(5) + y(1))) + sin(u(3))*(u(2) - (x(6) + y(2)))) + k(5)*(u(5) - (x(4) + y(5)));

%write control to car input
carInput(1) = delta;
carInput(2) = ax;

%------------- END OF CODE --------------