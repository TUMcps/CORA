function options = checkOptionsReach(obj,options)
% checkOptionsReach - checks if all necessary options
%   are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options)
%
% Inputs:
%    obj     - parallelHybridAutomaton object
%    options - options for parallelHybridAutomaton (1x1 struct)
%
% Outputs:
%    options - struct with adapted options
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      19-Feb-2019
% Last update:  04-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';

% check necessary options for reach for every component
numComp = length(obj.components);
for i=1:numComp
    % loop through every component
    comp = obj.components{i};
    whichCheck = check_flatHA_locations(options, comp);
end


% check options uniquely present in parallel hybrid automata...

% MANDATORY OPTIONS -------------------------------------------------------

check_R0(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
check_timeStep(options, obj,checkName);

options = check_pHA_inputs(options, obj);

check_guardIntersect(options, obj);
check_enclose(options, obj);
check_guardOrder(options, checkName);

options = check_startLoc_finalLoc(options, obj, 1);
options = check_reductionTechnique(options, obj);
options = check_linAlg(options, obj);

% OPTIONAL PARAMETERS -----------------------------------------------------

check_intersectInvariant(options, obj);


% WARNINGS ----------------------------------------------------------------

% extend validFields depending on contained obj (whichCheck)
validFields = getValidFields('parallelHybridAutomaton_reach');
if ~strcmp(options.guardIntersect, 'pancake') && ...
        ~strcmp(options.guardIntersect, 'hyperplaneMap')
    % enclose only if guardIntersect != pancake, hyperplaneMap
    validFields = [validFields; 'enclose'];
end

% extend validFields depending on systems in locations (whichCheck)
% validFieldsExt idx correspond to whichCheck

systems = {'linearSys'; 'linParamSys';...
    'nonlinDASys'; 'nonlinearSys'; 'nonlinParamSys'};
for i=1:length(whichCheck)
    if whichCheck(i)
        % concat nonredundant options for all contained systems
        validFieldsExt = getValidFields(systems{i});
        validFields = [validFields; validFieldsExt];
    end
end
    
printRedundancies(options, validFields);


end

%------------- END OF CODE --------------
