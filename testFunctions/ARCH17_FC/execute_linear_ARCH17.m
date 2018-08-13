function execute_linear_ARCH17()
% execute_linear_ARCH17 - executes all examples from the ARCH17 friendly 
% competition in the category linear continuous dynamics
%
% Syntax:  
%    execute_linear_ARCH17
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

disp('start building example')
example_linear_reach_ARCH17_building
disp('end building example')
disp('start platoon example: PLAD01/BND30')
example_linearParam_reach_ARCH17_platoon_deterministic_pass30
disp('end platoon example: PLAD01/BND30')
disp('start platoon example: PLAD01/BND42')
example_linearParam_reach_ARCH17_platoon_deterministic_pass42
disp('end platoon example: PLAD01/BND42')
disp('start platoon example: PLAA01/BND42')
example_linearParam_reach_ARCH17_platoon_pass42
disp('end platoon example: PLAA01/BND42')
disp('start platoon example: PLAA01/BND50')
example_linearParam_reach_ARCH17_platoon_pass50
disp('end platoon example: PLAA01/BND50')
disp('start platoon example: PLAN01/UND50')
example_linearParam_reach_ARCH17_platoon_unbounded
disp('end platoon example: PLAN01/UND50')
disp('start gearbox example')
example_linearHybrid_reach_ARCH17_gearbox
disp('end gearbox example')

%------------- END OF CODE --------------
