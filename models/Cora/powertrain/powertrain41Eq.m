function f = powertrain41Eq(x,u,p)
% powertrain41Eq - dynamic for the 41-dimensional power train system 
%                 (see Sec. 6 in [1])
%
% Syntax:  
%    f = powertrain41Eq(x,u,p)
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
f(6,1) = 1/p.J_l*(p.k_i*(x(40) - x(5)) - u(2) - p.b_l*x(6)); %\dot{Theta}_l
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
f(17,1) = p.J_i*(p.k_i*(x(14) - x(16)) - p.k_i*(x(16) - x(18)) - p.b_i*x(17)); %\dot{Theta}_5
f(18,1) = x(19); %Theta_6
f(19,1) = p.J_i*(p.k_i*(x(16) - x(18)) - p.k_i*(x(18) - x(20)) - p.b_i*x(19)); %\dot{Theta}_6
f(20,1) = x(21); %Theta_7
f(21,1) = p.J_i*(p.k_i*(x(18) - x(20)) - p.k_i*(x(20) - x(22)) - p.b_i*x(21)); %\dot{Theta}_7
f(22,1) = x(23); %Theta_8
f(23,1) = p.J_i*(p.k_i*(x(20) - x(22)) - p.k_i*(x(22) - x(24)) - p.b_i*x(23)); %\dot{Theta}_8
f(24,1) = x(25); %Theta_9
f(25,1) = p.J_i*(p.k_i*(x(22) - x(24)) - p.k_i*(x(24) - x(26)) - p.b_i*x(25)); %\dot{Theta}_9
f(26,1) = x(27); %Theta_10
f(27,1) = p.J_i*(p.k_i*(x(24) - x(26)) - p.k_i*(x(26) - x(28)) - p.b_i*x(27)); %\dot{Theta}_10
f(28,1) = x(29); %Theta_11
f(29,1) = p.J_i*(p.k_i*(x(26) - x(28)) - p.k_i*(x(28) - x(30)) - p.b_i*x(29)); %\dot{Theta}_11
f(30,1) = x(31); %Theta_12
f(31,1) = p.J_i*(p.k_i*(x(28) - x(30)) - p.k_i*(x(30) - x(32)) - p.b_i*x(31)); %\dot{Theta}_12
f(32,1) = x(33); %Theta_13
f(33,1) = p.J_i*(p.k_i*(x(30) - x(32)) - p.k_i*(x(32) - x(34)) - p.b_i*x(33)); %\dot{Theta}_13
f(34,1) = x(35); %Theta_14
f(35,1) = p.J_i*(p.k_i*(x(32) - x(34)) - p.k_i*(x(34) - x(36)) - p.b_i*x(35)); %\dot{Theta}_14
f(36,1) = x(37); %Theta_15
f(37,1) = p.J_i*(p.k_i*(x(34) - x(36)) - p.k_i*(x(36) - x(38)) - p.b_i*x(37)); %\dot{Theta}_15
f(38,1) = x(39); %Theta_16
f(39,1) = p.J_i*(p.k_i*(x(36) - x(38)) - p.k_i*(x(38) - x(40)) - p.b_i*x(39)); %\dot{Theta}_16
f(40,1) = x(41); %Theta_17
f(41,1) = p.J_i*(p.k_i*(x(38) - x(40)) - p.k_i*(x(40) - x(5)) - p.b_i*x(41)); %\dot{Theta}_17

%------------- END OF CODE --------------