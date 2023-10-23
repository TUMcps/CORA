function f = vehicleDynamics_KS_cog(x,u,p)
% vehicleDynamics_KS_cog - kinematic single-track vehicle dynamics 
% reference point: center of mass
%
% Inputs:
%    x - vehicle state vector
%    u - vehicle input vector
%    p - vehicle parameter structure
%
% Outputs:
%    f - right-hand side of differential equations
%
% Author:       Gerald Würsching
% Written:      17-November-2020
% Last update:  17-November-2020
% Last revision:---

%------------- BEGIN CODE --------------
%states
%x1 = s_x x-position in a global coordinate system
%x2 = s_y y-position in a global coordinate system
%x3 = δ steering angle of front wheels
%x4 = u velocity at center of mass
%x5 = Ψ yaw angle

%wheelbase
l_wb = p.a + p.b;

%consider steering constraints
u(1) = steeringConstraints(x(3),u(1),p.steering);

%consider acceleration constraints
u(2) = accelerationConstraints(x(4),u(2),p.longitudinal);

beta = atan(tan(x(3)) * p.b/l_wb) * sign(x(4));
%system dynamics
f(1,1) = x(4)*cos(beta + x(5));
f(2,1) = x(4)*sin(beta + x(5));
f(3,1) = u(1);
f(4,1) = u(2);
f(5,1) = x(4)*cos(beta)*tan(x(3))/l_wb;

%------------- END OF CODE --------------
