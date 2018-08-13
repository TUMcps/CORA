function execute_nonlinear_ARCH17()
% execute_nonlinear_ARCH17 - executes all examples from the ARCH17 friendly 
% competition in the category linear continuous dynamics
%
% Syntax:  
%    execute_nonlinear_ARCH17
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
% Written:      14-May-2017
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

disp('start van-der-Pol example')
example_nonlinear_reach_ARCH17_vanDerPol_pseudoInvariant
disp('end van-der-Pol example')
disp('start Laub-Loomis example with small initial set')
example_nonlinear_reach_ARCH17_sevenDim_small
disp('end Laub-Loomis example with small initial set')
disp('start Laub-Loomis example with large initial set')
example_nonlinear_reach_ARCH17_sevenDim
disp('end Laub-Loomis example with large initial set')
disp('start quadrocopter example')
example_nonlinear_reach_ARCH17_quadrocopterControlled
disp('end quadrocopter example')

%------------- END OF CODE --------------
