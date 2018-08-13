function execute_nonlinear_ARCH18()
% execute_nonlinear_ARCH18 - executes all examples from the ARCH18 friendly 
% competition in the category linear continuous dynamics
%
% Syntax:  
%    execute_nonlinear_ARCH18
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
% Author:       Niklas Kochdumper
% Written:      20-June-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

fprintf('\n\n\nBenchmark Van-der-Pol:\n\n')

disp('start van-der-Pol example')
example_nonlinear_reach_ARCH18_vanDerPol_zonoGirard
disp('end van-der-Pol example')



fprintf('\n\n\nBenchmark Laub-Loomis:\n\n')

disp('start Laub-Loomis example with small initial set')
example_nonlinear_reach_ARCH18_laubLoomis_small
disp('end Laub-Loomis example with small initial set')

disp('start Laub-Loomis example with large initial set')
example_nonlinear_reach_ARCH18_laubLoomis
disp('end Laub-Loomis example with large initial set')



fprintf('\n\n\nBenchmark Quadrocopter:\n\n')

disp('start quadrocopter example')
example_nonlinear_reach_ARCH18_quadrocopterControlled
disp('end quadrocopter example')



fprintf('\n\n\nBenchmark Spacecraft Rendezvous:\n\n')

disp('start spacecraft example')
example_nonlinear_reach_ARCH18_rendezvousSX4p_zonoGirard
disp('end spacecraft example')

%------------- END OF CODE --------------