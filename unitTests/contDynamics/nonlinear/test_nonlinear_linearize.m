function res = test_nonlinear_linearize(~)
% test_nonlinear_linearize - unit_test_function of linearizing nonlinear 
% dynamics
%
% Checks the linearization of the nonlinearSys class for the 6 tank example;
% It is checked whether the A and B matrix are correct for a particular
% linearization point
%
% Syntax:  
%    res = test_nonlinear_linearize(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      30-July-2017
% Last update:  12-September-2017
% Last revision:---


%------------- BEGIN CODE --------------

dim=6;

%set options --------------------------------------------------------------
options.uTrans = 0;
options.U = zonotope([0,0.005]); %input for reachability analysis
options.tensorOrder = 2;
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
tank = nonlinearSys(6,1,@tank6Eq,options); %initialize tank system
%--------------------------------------------------------------------------

%linearize system
R0=zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim)]); %initial state for reachability analysis
[~,linSys,linOptions] = linearize(tank,options,R0);

% provide ground truth-----------------------------------------------------
A_true = [  -0.023490689645049, 0, 0, 0, 0, -0.010000000000000; ...
            0.023490689645049, -0.016610425942763, 0, 0, 0, 0; ...
            0, 0.016610425942763, -0.016610425942763, 0, 0, 0; ...
            0, 0, 0.016610425942763, -0.023490689645049, 0, 0; ...
            0, 0, 0, 0.023490689645049, -0.010505355776936, 0; ...
            0, 0, 0, 0, 0.010505355776936, -0.016610425942763];
U_true_center = zeros(dim,1);
U_true_generator = [0.005; zeros(dim-1,1)];
%--------------------------------------------------------------------------

%compare with obtained values
res_1 = (max(max(abs(linSys.A - A_true))) <= 1e-12);
res_2 = (max(abs(center(linOptions.U) - U_true_center)) <= 1e-12);
res_3 = (max(abs(generators(linOptions.U) - U_true_generator)) <= 1e-12);

%final result
res = res_1*res_2*res_3;

%------------- END OF CODE --------------
