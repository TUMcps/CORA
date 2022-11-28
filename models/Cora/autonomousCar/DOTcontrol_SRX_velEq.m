function [carInput,k] = DOTcontrol_SRX_velEq(x,uComb)
% DOTcontrol_SRX_vel - provides the steering angle speed of the Cadillac 
%                      SRX (see [1])
%
% Syntax:  
%    [carInput,k] = DOTcontrol_SRX_velEq(x,uComb)
%
% Inputs:
%    x - state
%    uComb - combination of reference trajectory u_ref and noise y 
%
% Outputs:
%    carInput - steering angle for the car
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
% Written:      01-March-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% %control parameters
% k(1) = 0.2;
% k(2) = 2;
% k(3) = 0.3;
% k(6) = 1; %gain for steering wheel
%control parameters
k(1) = 0.5;
k(2) = 3;
k(3) = 1;
k(6) = 2; %gain for steering wheel
%k(6) = 1/3.4*2; %gain for steering wheel

% separate in reference u and noise y
u(1:5) = uComb(1:5);
y(1:4) = uComb(5+(1:4));

%compute steering angle and acceleration input from reference trajectory
% steering angle delta
delta = k(1)*(-sin(u(3))*(u(1) - (x(4) + y(1))) + cos(u(3))*(u(2) - (x(5) + y(2)))) ...
    + k(2)*(u(3) - (x(2) + y(3))) + k(3)*(u(4) - (x(3) + y(4)));
%steering wheel speed
delta_dot = k(6)*(delta - x(6));


%write control to car input
carInput(1) = delta_dot;
carInput(2) = u(5);

%------------- END OF CODE --------------