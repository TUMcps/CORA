function options = checkOptionsSimRandom(obj,options)
% checkOptionsSimRandom - checks if all necessary options
%                         are there and have valid values
%
% Syntax:
%   option = checkOptionsSimRandom(opt)
%
% Inputs:
%    obj     - parallelHybridAutomaton object
%    options - options for parallelHybridAutomaton simulation (1x1 struct)
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

checkName = "checkOptionsSimRandom";


% MANDATORY OPTIONS -------------------------------------------------------

check_metasim(options,obj,true);

check_R0(options, obj);
options = check_pHA_inputs(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
options = check_startLoc_finalLoc(options, obj, 1);


% WARNINGS ----------------------------------------------------------------

% necessary fields 
validFields = getValidFields('parallelHybridAutomaton_simRand');
printRedundancies(options, validFields);

end

%------------- END OF CODE --------------
