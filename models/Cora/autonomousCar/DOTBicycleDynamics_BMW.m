function f = DOTBicycleDynamics_BMW(x,u)
% DOTBicycleDynamics_BMW - generates bicycle model for a BMW 
%
% Syntax:  
%    f = DOTBicycleDynamics_BMW(t,x,u)
%
% Inputs:
%    x - state vector
%    u - input vector (here: reference trajectory)
%
% Outputs:
%    f - time-derivative of the state vector
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      26-August-2011
% Last update:  30-August-2011
%               14-June-2023
% Last revision:---

%------------- BEGIN CODE --------------

%load parameters
g = 9.81; %[m/s^2]

%get model paramters
p = BMWparameters();

%create equivalent bicycle parameters
mu = p.mu;
C_Sf = -p.tire.p_ky1/p.tire.p_dy1; 
C_Sr = -p.tire.p_ky1/p.tire.p_dy1; 
lf = p.a;
lr = p.b;
h = p.h_s;
m = p.m;
I = p.I_z;

%states
%x1 = β slip angle at vehicle center
%x2 = Ψ yaw angle
%x3 = Ψ yaw rate
%x4 = u velocity in x-direction
%x5 = s_x x-position in a global coordinate system
%x6 = s_y y-position in a global coordinate system

%u1 = delta_w steering angle of front wheels
%u2 = ax longitudinal acceleration


%system dynamics
f(1,1) = (mu/(x(4)^2*(lr+lf))*(C_Sr*(g*lf + u(2)*h)*lr - C_Sf*(g*lr - u(2)*h)*lf)-1)*x(3) ...
    -mu/(x(4)*(lr+lf))*(C_Sr*(g*lf + u(2)*h) + C_Sf*(g*lr-u(2)*h))*x(1) ...
    +mu/(x(4)*(lr+lf))*(C_Sf*(g*lr-u(2)*h))*u(1);
f(2,1) = x(3);
f(3,1) = -mu*m/(x(4)*I*(lr+lf))*(lf^2*C_Sf*(g*lr-u(2)*h) + lr^2*C_Sr*(g*lf + u(2)*h))*x(3) ...
    +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf + u(2)*h) - lf*C_Sf*(g*lr - u(2)*h))*x(1) ...
    +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr - u(2)*h)*u(1);
f(4,1) = u(2);
f(5,1) = x(4)*cos(x(1) + x(2));
f(6,1) = x(4)*sin(x(1) + x(2));


%------------- END OF CODE --------------