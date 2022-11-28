function f = powertrain101Eq(x,u,p)
% powertrain101Eq - dynamic for the 101-dimensional power train system 
%                   (see Sec. 6 in [1])
%
% Syntax:  
%    f = powertrain101Eq(x,u,p)
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
f(6,1) = 1/p.J_l*(p.k_i*(x(100) - x(5)) - u(2) - p.b_l*x(6)); %\dot{Theta}_l
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
f(41,1) = p.J_i*(p.k_i*(x(38) - x(40)) - p.k_i*(x(40) - x(42)) - p.b_i*x(41)); %\dot{Theta}_17
f(42,1) = x(43); %Theta_18
f(43,1) = p.J_i*(p.k_i*(x(40) - x(42)) - p.k_i*(x(42) - x(44)) - p.b_i*x(43)); %\dot{Theta}_18
f(44,1) = x(45); %Theta_19
f(45,1) = p.J_i*(p.k_i*(x(42) - x(44)) - p.k_i*(x(44) - x(46)) - p.b_i*x(45)); %\dot{Theta}_19
f(46,1) = x(47); %Theta_20
f(47,1) = p.J_i*(p.k_i*(x(44) - x(46)) - p.k_i*(x(46) - x(48)) - p.b_i*x(47)); %\dot{Theta}_20
f(48,1) = x(49); %Theta_21
f(49,1) = p.J_i*(p.k_i*(x(46) - x(48)) - p.k_i*(x(48) - x(50)) - p.b_i*x(49)); %\dot{Theta}_21
f(50,1) = x(51); %Theta_22
f(51,1) = p.J_i*(p.k_i*(x(48) - x(50)) - p.k_i*(x(50) - x(52)) - p.b_i*x(51)); %\dot{Theta}_22
f(52,1) = x(53); %Theta_23
f(53,1) = p.J_i*(p.k_i*(x(50) - x(52)) - p.k_i*(x(52) - x(54)) - p.b_i*x(53)); %\dot{Theta}_23
f(54,1) = x(55); %Theta_24
f(55,1) = p.J_i*(p.k_i*(x(52) - x(54)) - p.k_i*(x(54) - x(56)) - p.b_i*x(55)); %\dot{Theta}_24
f(56,1) = x(57); %Theta_25
f(57,1) = p.J_i*(p.k_i*(x(54) - x(56)) - p.k_i*(x(56) - x(58)) - p.b_i*x(57)); %\dot{Theta}_25
f(58,1) = x(59); %Theta_26
f(59,1) = p.J_i*(p.k_i*(x(56) - x(58)) - p.k_i*(x(58) - x(60)) - p.b_i*x(59)); %\dot{Theta}_26
f(60,1) = x(61); %Theta_27
f(61,1) = p.J_i*(p.k_i*(x(58) - x(60)) - p.k_i*(x(60) - x(62)) - p.b_i*x(61)); %\dot{Theta}_27
f(62,1) = x(63); %Theta_28
f(63,1) = p.J_i*(p.k_i*(x(60) - x(62)) - p.k_i*(x(62) - x(64)) - p.b_i*x(63)); %\dot{Theta}_28
f(64,1) = x(65); %Theta_29
f(65,1) = p.J_i*(p.k_i*(x(62) - x(64)) - p.k_i*(x(64) - x(66)) - p.b_i*x(65)); %\dot{Theta}_29
f(66,1) = x(67); %Theta_30
f(67,1) = p.J_i*(p.k_i*(x(64) - x(66)) - p.k_i*(x(66) - x(68)) - p.b_i*x(67)); %\dot{Theta}_30
f(68,1) = x(69); %Theta_31
f(69,1) = p.J_i*(p.k_i*(x(66) - x(68)) - p.k_i*(x(68) - x(70)) - p.b_i*x(69)); %\dot{Theta}_31
f(70,1) = x(71); %Theta_32
f(71,1) = p.J_i*(p.k_i*(x(68) - x(70)) - p.k_i*(x(70) - x(72)) - p.b_i*x(71)); %\dot{Theta}_32
f(72,1) = x(73); %Theta_33
f(73,1) = p.J_i*(p.k_i*(x(70) - x(72)) - p.k_i*(x(72) - x(74)) - p.b_i*x(73)); %\dot{Theta}_33
f(74,1) = x(75); %Theta_34
f(75,1) = p.J_i*(p.k_i*(x(72) - x(74)) - p.k_i*(x(74) - x(76)) - p.b_i*x(75)); %\dot{Theta}_34
f(76,1) = x(77); %Theta_35
f(77,1) = p.J_i*(p.k_i*(x(74) - x(76)) - p.k_i*(x(76) - x(78)) - p.b_i*x(77)); %\dot{Theta}_35
f(78,1) = x(79); %Theta_36
f(79,1) = p.J_i*(p.k_i*(x(76) - x(78)) - p.k_i*(x(78) - x(80)) - p.b_i*x(79)); %\dot{Theta}_36
f(80,1) = x(81); %Theta_37
f(81,1) = p.J_i*(p.k_i*(x(78) - x(80)) - p.k_i*(x(80) - x(82)) - p.b_i*x(81)); %\dot{Theta}_37
f(82,1) = x(83); %Theta_38
f(83,1) = p.J_i*(p.k_i*(x(80) - x(82)) - p.k_i*(x(82) - x(84)) - p.b_i*x(83)); %\dot{Theta}_38
f(84,1) = x(85); %Theta_39
f(85,1) = p.J_i*(p.k_i*(x(82) - x(84)) - p.k_i*(x(84) - x(86)) - p.b_i*x(85)); %\dot{Theta}_39
f(86,1) = x(87); %Theta_40
f(87,1) = p.J_i*(p.k_i*(x(84) - x(86)) - p.k_i*(x(86) - x(88)) - p.b_i*x(87)); %\dot{Theta}_40
f(88,1) = x(89); %Theta_41
f(89,1) = p.J_i*(p.k_i*(x(86) - x(88)) - p.k_i*(x(88) - x(90)) - p.b_i*x(89)); %\dot{Theta}_41
f(90,1) = x(91); %Theta_42
f(91,1) = p.J_i*(p.k_i*(x(88) - x(90)) - p.k_i*(x(90) - x(92)) - p.b_i*x(91)); %\dot{Theta}_42
f(92,1) = x(93); %Theta_43
f(93,1) = p.J_i*(p.k_i*(x(90) - x(92)) - p.k_i*(x(92) - x(94)) - p.b_i*x(93)); %\dot{Theta}_43
f(94,1) = x(95); %Theta_44
f(95,1) = p.J_i*(p.k_i*(x(92) - x(94)) - p.k_i*(x(94) - x(96)) - p.b_i*x(95)); %\dot{Theta}_44
f(96,1) = x(97); %Theta_45
f(97,1) = p.J_i*(p.k_i*(x(94) - x(96)) - p.k_i*(x(96) - x(98)) - p.b_i*x(97)); %\dot{Theta}_45
f(98,1) = x(99); %Theta_46
f(99,1) = p.J_i*(p.k_i*(x(96) - x(98)) - p.k_i*(x(98) - x(100)) - p.b_i*x(99)); %\dot{Theta}_46
f(100,1) = x(101); %Theta_47
f(101,1) = p.J_i*(p.k_i*(x(98) - x(100)) - p.k_i*(x(100) - x(5)) - p.b_i*x(101)); %\dot{Theta}_47

%------------- END OF CODE --------------