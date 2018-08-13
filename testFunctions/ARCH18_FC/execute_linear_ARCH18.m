function execute_linear_ARCH18()
% execute_linear_ARCH18 - executes all examples from the ARCH18 friendly 
% competition in the category linear continuous dynamics
%
% Syntax:  
%    execute_linear_ARCH18
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
% Written:      20-April-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

fprintf('\n\n\nBenchmark ISS:\n\n')

disp('start ISS example: ISSF01/ISS01')
example_linear_reach_ARCH18_iss_safe
disp('end ISS example: ISSF01/ISS01')

disp('start ISS example: ISSF01/ISU01')
example_linear_reach_ARCH18_iss_unsafe
disp('end ISS example: ISSF01/ISU01')

disp('start ISS example: ISSC01/ISS02')
example_linear_reach_ARCH18_iss_Uconst_safe
disp('end ISS example: ISSC01/ISS02')

disp('start ISS example: ISSC01/ISU02')
example_linear_reach_ARCH18_iss_Uconst_unsafe
disp('end ISS example: ISSC01/ISU02')



fprintf('\n\n\nBenchmark Spacecraft Rendezvous:\n\n')

disp('start spacecraft example: SRNA01/SR01')
example_linear_reach_ARCH18_rendezvousSX4np_zonoGirard
disp('end spacecraft example: SRNA01/SR01')

disp('start spacecraft example: SRA01/SR02')
example_linear_reach_ARCH18_rendezvousSX4p_zonoGirard
disp('end spacecraft example: SRA01/SR02')



fprintf('\n\n\nBenchmark Drivetrain:\n\n')

disp('start drivetrain example: DTN01')
example_hybrid_reach_ARCH18_powerTrain_DTN01
disp('end drivetrain example: DTN01')

disp('start drivetrain example: DTN02')
example_hybrid_reach_ARCH18_powerTrain_DTN02
disp('end drivetrain example: DTN02')

disp('start drivetrain example: DTN03')
example_hybrid_reach_ARCH18_powerTrain_DTN03
disp('end drivetrain example: DTN03')

disp('start drivetrain example: DTN04')
example_hybrid_reach_ARCH18_powerTrain_DTN04
disp('end drivetrain example: DTN04')

disp('start drivetrain example: DTN05')
example_hybrid_reach_ARCH18_powerTrain_DTN05
disp('end drivetrain example: DTN05')

disp('start drivetrain example: DTN06')
example_hybrid_reach_ARCH18_powerTrain_DTN06
disp('end drivetrain example: DTN06')



fprintf('\n\n\nBenchmark Building:\n\n')

disp('start building example: BLDC01/BDS01')
example_linear_reach_ARCH18_building_constInput
disp('end building example: BLDC01/BDS01')

disp('start building example: BLDF01/BDS01')
example_linear_reach_ARCH18_building
disp('end building example: BLDF01/BDS01')



fprintf('\n\n\nBenchmark Platooning:\n\n')

disp('start platoon example: PLAD01/BND30')
example_linearParam_reach_ARCH18_platoon_deterministic_pass30
disp('end platoon example: PLAD01/BND30')

disp('start platoon example: PLAD01/BND42')
example_linearParam_reach_ARCH18_platoon_deterministic_pass42
disp('end platoon example: PLAD01/BND42')

disp('start platoon example: PLAA01/BND42')
example_linearParam_reach_ARCH18_platoon_pass42
disp('end platoon example: PLAA01/BND42')

disp('start platoon example: PLAA01/BND50')
example_linearParam_reach_ARCH18_platoon_pass50
disp('end platoon example: PLAA01/BND50')

disp('start platoon example: PLAN01/UND50')
example_linearParam_reach_ARCH18_platoon_unbounded
disp('end platoon example: PLAN01/UND50')



fprintf('\n\n\nBenchmark Gearbox:\n\n')

disp('start gearbox example')
example_linearHybrid_reach_ARCH18_gearbox_zonoGirard
disp('end gearbox example')

%------------- END OF CODE --------------
