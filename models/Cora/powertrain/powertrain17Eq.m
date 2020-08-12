function f = powertrain17Eq(x,u,p)
% powertrain17Eq - dynamic for the 17-dimensional power train system 
%                 (see Sec. 6 in [1])
%
% Syntax:  
%    f = powertrain17Eq(x,u,p)
%
% Inputs:
%    x - state vector
%    u - input vector
%    p - struct storing the model parameter
%
% Outputs:
%    f - time-derivate of the system state
%
% References:
%   [1] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems"

% Author:       Matthias Althoff
% Written:      21-September-2011
% Last update:  23-December-2019
% Last revision:---

%------------- BEGIN CODE --------------

%control
v = p.k_K*(p.i*x(4) - x(7)) ...
    + p.k_KD*(p.i*u(1) - 1/p.J_m*(x(2) - 1/p.i*p.k*(x(1) - p.alpha) - p.b_m*x(7))) ...
    + p.k_KI*(p.i*x(3) - p.i*(x(1) + x(8)));


%plant model
f(1,1) = 1/p.i*x(7) - x(9); %Theta_d
f(2,1) = (v - x(2))/p.tau_eng; %T_m
f(3,1) = x(4); %Theta_ref
f(4,1) = u(1); %\dot{Theta}_ref
f(5,1) = x(6); %Theta_l
f(6,1) = 1/p.J_l*(p.k_i*(x(16) - x(5)) - u(2) - p.b_l*x(6)); %\dot{Theta}_l
f(7,1) = 1/p.J_m*(x(2) - 1/p.i*p.k*(x(1) - p.alpha) - p.b_m*x(7)); %\dot{Theta}_m
f(8,1) = x(9); %Theta_1
f(9,1) = p.J_i*(p.k*(x(1) - p.alpha) - p.k_i*(x(8) - x(10)) - p.b_i*x(9)); %\dot{Theta}_1
f(10,1) = x(11); %Theta_2
f(11,1) = p.J_i*(p.k_i*(x(8) - x(10)) - p.k_i*(x(10) - x(12)) - p.b_i*x(11)); %\dot{Theta}_2
f(12,1) = x(13); %Theta_3
f(13,1) = p.J_i*(p.k_i*(x(10) - x(12)) - p.k_i*(x(12) - x(14)) - p.b_i*x(13)); %\dot{Theta}_3
f(14,1) = x(15); %Theta_4
f(15,1) = p.J_i*(p.k_i*(x(12) - x(14)) - p.k_i*(x(14) - x(16)) - p.b_i*x(15)); %\dot{Theta}_4
f(16,1) = x(17); %Theta_5
f(17,1) = p.J_i*(p.k_i*(x(14) - x(16)) - p.k_i*(x(16) - x(5)) - p.b_i*x(17)); %\dot{Theta}_5

%------------- END OF CODE --------------