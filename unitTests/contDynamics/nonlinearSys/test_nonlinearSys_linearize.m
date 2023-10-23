function res = test_nonlinearSys_linearize
% test_nonlinearSys_linearize - unit_test_function of linearizing nonlinear 
%    dynamics: Checks the linearization of the nonlinearSys class
%    for the 6 tank example; It is checked whether the A and B matrix
%    are correct for a particular linearization point
%
% Syntax:
%    res = test_nonlinearSys_linearize
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       30-July-2017
% Last update:   12-September-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dim_x=6;

% Reachability Settings ---------------------------------------------------

options.U = zonotope([0,0.005]); %input for reachability analysis
options.uTrans = 0; %since only linearize is called
options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------

tank = nonlinearSys(@tank6Eq); %initialize tank system


% linearize system
R0 = zonotope([[2; 4; 4; 2; 10; 4],0.2*eye(dim_x)]); %initial state for reachability analysis
[~,linSys,linOptions] = linearize(tank,options,R0);

% provide ground truth ----------------------------------------------------
A_true = [  -0.023490689645049, 0, 0, 0, 0, -0.010000000000000; ...
            0.023490689645049, -0.016610425942763, 0, 0, 0, 0; ...
            0, 0.016610425942763, -0.016610425942763, 0, 0, 0; ...
            0, 0, 0.016610425942763, -0.023490689645049, 0, 0; ...
            0, 0, 0, 0.023490689645049, -0.010505355776936, 0; ...
            0, 0, 0, 0, 0.010505355776936, -0.016610425942763];
U_true_center = zeros(dim_x,1);
U_true_generator = [0.005; zeros(dim_x-1,1)];
% -------------------------------------------------------------------------

%compare with obtained values
res_1 = (max(max(abs(linSys.A - A_true))) <= 1e-12);
res_2 = (max(abs(center(linOptions.U) - U_true_center)) <= 1e-12);
res_3 = (max(abs(generators(linOptions.U) - U_true_generator)) <= 1e-12);

%final result
res = res_1 && res_2 && res_3;

% ------------------------------ END OF CODE ------------------------------
