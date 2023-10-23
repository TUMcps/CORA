function powerSystem2cora(filename, varargin)
% powerSystem2cora - convert a power system description to a CORA model
%
% Syntax:
%    powerSystem2cora(filename);
%
% Inputs:
%    filename - name of the CommonRoad XML-file scenario (specified as string)
%    'verbose', true/false - name-value pair specifying whether conversion
%           information should be printed to the terminal. Default is true.
%
%    structure:
%    .name - name of power system
%    .Y - admittance matrix
%    .Pd - active power demand
%    .Qd - reactive power demand
%    .VM - voltage magnitude
%    .genParam - generator parameters
%    .bus.output - vector of output buses
%    .bus.input - vector of input buses
%    .bus.generator - vector of generator buses
%    .bus.load - vector of load buses
%    .bus.fault - index of faulty bus
%    .bus.slack - index of slack bus
%
% Outputs:
%    x - ...
%    y - ...
%
% Example:
%    [x,y] = powerSystem2cora('IEEE30');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: spaceex2cora

% Authors:       Matthias Althoff
% Written:       14-April-2022
% Last update:   25-April-2023 (MW, change destination folder)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load scenario
scenario = load(filename);

% reorder bus indices (generatorBuses, loadBuses, inputBuses)
scenario = reorderBuses(scenario); 

% create power system variables
[powVariables, NrOf] = symPowerVariables(scenario);

% constraint equations
disp('create constraint equations');
generateConstraintEquations(scenario,powVariables);

% dynamic equations
disp('create dynamic equations');
generateDynamicEquations(scenario,powVariables);
    
    
% create CORA DAE model
% name of model
name = [scenario.name,'_model'];
% create string
str = [name,' = nonlinDASys(''',scenario.name,''',@',[scenario.name,'_dyn'],...
    ',@',[scenario.name,'_con'],',',num2str(NrOf.states),',',num2str(NrOf.inputs),',',num2str(NrOf.constraints),');'];
%evaluate string
eval(str);

%% save model
% set path
path = [CORAROOT filesep 'models' filesep 'powerSystemsConverted'];
if ~isfolder(path)
    mkdir(path);
end
save([path filesep name], name);
% remove and add path so that file can be found
warOrig = warning;
warning('off','all');
rmpath(path);
warning(warOrig);
addpath(path);

% ------------------------------ END OF CODE ------------------------------
