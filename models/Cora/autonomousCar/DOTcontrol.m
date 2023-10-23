function [carInput,k] = DOTcontrol(x,varargin)
% DOTcontrol - provides the steering angle and acceleration of the vehicle
%              (adapted from [1] and merged with vehicle model in CommonRoad)
%
% Syntax:  
%    f = DOTcontrol(x,u,y)
%
% Inputs:
%    x - state vector
%    uComb - combination of reference trajectory u_ref and noise y 
%
% Outputs:
%    carInput - control inputs for the car
%    k - control parameter
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DOTBicycleDynamics_controlled_BMW

% Author:       Matthias Althoff
% Written:      23-August-2011
% Last update:  15-June-2023
%               22-June-2023
%               21-September-2023 (TL, parse input)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input
if nargin < 2
    throw(CORAerror("CORA:notEnoughInputArgs", 2))
elseif nargin == 2
    uComb = varargin{1};
elseif nargin == 3
    uComb = [varargin{1};varargin{2}];
else
    throw(CORAerror("CORA:tooManyInputArgs", 3))
end


%control parameters
k(1) = 0.2;
k(2) = 2;
k(3) = 0.3;
k(4) = 1;
k(5) = 10; 
k(6) = 10;


%states
%x1 = s_x x-position in a global coordinate system
%x2 = s_y y-position in a global coordinate system
%x3 = δ steering angle of front wheels
%x4 = u velocity in x-direction
%x5 = Ψ yaw angle
%x6 = Ψ yaw rate
%x7 = β slip angle at vehicle center

% controller outputs
% carInput1 = delta_dot steering velocity of front wheels
% carInput2 = ax longitudinal acceleration

% separate in reference u and noise y
% u1: sx_d desired x-position in a global coordinate system
% u2: sy_d desired y-position in a global coordinate system
% u3: Psi_d desired yaw angle
% u4: \dot{Psi}_d desired yaw rate
% u5: v desired velocity

u(1:5) = uComb(1:5); % reference inputs
v(1:5) = uComb(5+(1:5)); % sensor noise v

%compute steering angle and acceleration input from reference trajectory
% steering angle delta
delta_d = k(1)*(-sin(u(3))*(u(1) - (x(1) + v(1))) + cos(u(3))*(u(2) - (x(2) + v(2)))) ...
    + k(2)*(u(3) - (x(5) + v(3))) + k(3)*(u(4) - (x(6) + v(4))); % desired steeing angle
delta_dot = k(6)*(delta_d - x(3)); % P controller to obtaindesired steering angle
% longitudinal acceleration ax 
ax = k(4)*(cos(u(3))*(u(1) - (x(1) + v(1))) + sin(u(3))*(u(2) - (x(2) + v(2)))) + k(5)*(u(5) - (x(4) + v(5)));

%write control to car input
carInput(1,1) = delta_dot;
carInput(2,1) = ax;

%------------- END OF CODE --------------