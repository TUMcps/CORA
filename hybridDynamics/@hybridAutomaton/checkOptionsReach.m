function options = checkOptionsReach(obj,options)
% checkOptionsReach - checks if all necessary options
%   are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options)
%
% Inputs:
%    obj     - hybridAutomaton object (unused)
%    options - options for hybridAutomaton (1x1 struct)
%
% Outputs:
%    options - options for hybridAutomaton (1x1 struct)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      10-Feb-2019
% Last update:  09-Sep-2019
%               04-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';

% information for validFields of redundancy warnings
whichCheck = check_flatHA_locations(options, obj);

% check options uniquely present in hybrid automata...

% MANDATORY OPTIONS -------------------------------------------------------

check_flatHA_output(obj);

options = check_linAlg(options, obj);
options = check_flatHA_inputs(options, obj);
options = check_timeStepLoc(options, obj);

check_guardIntersect(options, obj);
check_enclose(options, obj);
check_guardOrder(options, checkName);

options = check_startLoc_finalLoc(options, obj, 0);


% OPTIONAL PARAMETERS -----------------------------------------------------

check_intersectInvariant(options, checkName);
check_compTimePoint(options, checkName);


% WARNINGS ----------------------------------------------------------------

validFields = getValidFields('hybridAutomaton_reach');
if ~strcmp(options.guardIntersect, 'pancake') && ...
    ~strcmp(options.guardIntersect, 'hyperplaneMap') && ...
    ~strcmp(options.guardIntersect,'levelSet')
    % enclose only if guardIntersect != pancake, hyperplaneMap, levelSet
    validFields = [validFields; 'enclose'];
end
if ismember(options.guardIntersect,{'hyperplaneMap','conZonotope'})
    validFields = [validFields; 'guardOrder'];
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

unique(validFields);
printRedundancies(options, validFields);


% default options
if ~isfield(options,'tStart')
    options.tStart = 0;
end

if ~isfield(options,'reductionTechnique')
    options.reductionTechnique = 'girard'; 
end

if any(whichCheck(3:end))
   if ~isfield(options,'maxError')
      options.maxError = Inf*ones(dim(options.R0),1); 
   end
   if ~isfield(options,'reductionInterval')
      options.reductionInterval = Inf; 
   end
end

if any(whichCheck(4:end))
    options = check_polyZono(options, obj);
end

end

%------------- END OF CODE --------------
