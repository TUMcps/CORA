function options = checkOptionsSimulate(obj,options)
% checkOptionsSimulate - checks if all necessary options
%   are there and have valid values
%
% Syntax:  
%    options = checkOptionsSimulate(obj,options)
%
% Inputs:
%    obj     - hybridAutomaton object
%    options - options for hybridAutomaton simulation (1x1 struct)
%
% Outputs:
%    options - options for hybridAutomaton simulation (1x1 struct)
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
% Written:      19-Feb-2019
% Last update:  09-Sep-2019
% Last revision:---

%------------- BEGIN CODE --------------

checkName = "checkOptionsSimulate";


% MANDATORY OPTIONS -------------------------------------------------------

locations = obj.location;
startLocation = locations{options.startLoc};
locationSys = startLocation.contDynamics;

% x0 has to exist, size has to match dim of startLoc
if ~isfield(options,'x0')
    printOptionMissing(checkName, 'x0');
elseif length(options.x0) ~= locationSys.dim
    printOptionOutOfRange(checkName, 'x0','params');
end


options = check_flatHA_inputsSim(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
options = check_startLoc_finalLoc(options, obj, 0);


% WARNINGS ----------------------------------------------------------------

% necessary fields 
validFields = getValidFields('hybridAutomaton_sim');
printRedundancies(options, validFields);


end

%------------- END OF CODE --------------
