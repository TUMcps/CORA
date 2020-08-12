function options = checkOptionsSimulate(obj,options)
% checkOptionsSimulate - checks if all necessary options
%   are there and have valid values
%
% Syntax:
%    options = checkOptionsSimulate(opt)
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

checkName = "checkOptionsSimulate";


% MANDATORY OPTIONS -------------------------------------------------------

check_pha_x0(options, obj);
options = check_pHA_inputsSim(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
options = check_startLoc_finalLoc(options, obj, 1);


% WARNINGS ----------------------------------------------------------------

% necessary fields 
validFields = getValidFields('parallelHybridAutomaton_sim');
printRedundancies(options, validFields);


end

%------------- END OF CODE --------------
