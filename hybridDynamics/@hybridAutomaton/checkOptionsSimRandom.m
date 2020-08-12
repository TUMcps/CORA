function options = checkOptionsSimRandom(obj,options)
% checkOptionsSimRandom - checks if all necessary options
%                         are there and have valid values
%
% Syntax:  
%    options = checkOptionsSimRandom(obj,options)
%
% Inputs:
%    obj     - hybridAutomaton object (unused)
%    options - options for hybridAutomaton random simulation (1x1 struct)
%
% Outputs:
%    options - options for hybridAutomaton random simulation (1x1 struct)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Niklas Kochdumper
% Written:      23-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

checkName = "checkOptionsReach";


% MANDATORY OPTIONS -------------------------------------------------------

check_metasim(options, obj, true);

locations = obj.location;
startLocation = locations{options.startLoc};
locationSys = startLocation.contDynamics;

% R0 has to exist, size has to match dim of startLoc
if ~isfield(options,'R0')
    printOptionMissing(obj, 'R0');
elseif dim(options.R0) ~= locationSys.dim
    printOptionOutOfRange(obj, 'R0');
end

options = check_flatHA_inputs(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
options = check_startLoc_finalLoc(options, obj, 0);


% WARNINGS ----------------------------------------------------------------

% necessary fields 
validFields = getValidFields('hybridAutomaton_simRand');
printRedundancies(options, validFields);


end

%------------- END OF CODE --------------
