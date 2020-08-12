function f = DOTBicycleDynamics_SRX_velEq(x,u)
% DOTBicycleDynamics_SRX_vel - generates bicycle model for the Cadillac SRX
%                              for a given velocity (see Eq. (1) in [1])
%
% Syntax:  
%    f = DOTBicycleDynamics_SRX_vel(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivative of the state vector
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
% See also: example_nonlinear_reach05_autonomousCar

% Author:       Matthias Althoff
% Written:      01-March-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load parameters
g = 9.81; %[m/s^2]

%get model paramters
p = SRXparameters();

%create equivalent bicycle parameters
%mu = p.tire.p_dy1;
mu = 1; %<-- hard coded
C_Sf = -p.tire.p_ky1/p.tire.p_dy1; %<--corrected
C_Sr = -p.tire.p_ky1/p.tire.p_dy1; %<--corrected
lf = p.a;
lr = p.b;
m = p.m;
I = p.I_z;
k_P = 1;
k_D = 1;
I_wheel = 3;
T_steer = 1;

%states
%x1 = β slip angle at vehicle center
%x2 = Ψ yaw angle
%x3 = Ψ yaw rate
%x5 = s_x x-position in a global coordinate system
%x6 = s_y y-position in a global coordinate system
%x7 = delta steering angle of front wheels

%u1 = delta_w steering angle velocity of front wheels
%u2 = velocity


%system dynamics
f(1,1) = (mu/(u(2)^2*(lr+lf))*(C_Sr*(g*lf)*lr - C_Sf*(g*lr)*lf)-1)*x(3) ...
    -mu/(u(2)*(lr+lf))*(C_Sr*(g*lf) + C_Sf*(g*lr))*x(1) ...
    +mu/(u(2)*(lr+lf))*(C_Sf*(g*lr))*x(6);
f(2,1) = x(3);
f(3,1) = -mu*m/(u(2)*I*(lr+lf))*(lf^2*C_Sf*(g*lr) + lr^2*C_Sr*(g*lf))*x(3) ...
    +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf) - lf*C_Sf*(g*lr))*x(1) ...
    +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr)*x(6);
f(4,1) = u(2)*cos(x(1) + x(2));
f(5,1) = u(2)*sin(x(1) + x(2));
f(6,1) = u(1);
%f(6,1) = -1/T_steer*x(6) + u(1);
% f(6,1) = x(7) + k_D/I_wheel*(u(1) - x(7));
% f(7,1) = k_P/I_wheel*(u(1) - x(7));


%------------- END OF CODE --------------