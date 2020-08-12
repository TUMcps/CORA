function options = checkOptionsReach(obj,options,hyb)
% checkOptionsReach - checks if all necessary options
%    are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options,hyb)
%
% Inputs:
%    obj     - linParamSys object (unused)
%    options - options for linearSys (1x1 struct)
%    hyb     - (optional) called from hybrid Automaton (0/1)
%
% Outputs:
%    options - options for linearSys (1x1 struct)
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
% Written:      17-Feb-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';
if nargin ~= 3
    hyb = 0;
end
suppressWarnings = hyb;

% DEFAULT OPTIONS ---------------------------------------------------------

options = check_reductionTechnique(options, obj);
options = check_verbose(options, obj);
options = check_compTimePoint(options, obj);

% MANDATORY OPTIONS -------------------------------------------------------

check_R0(options, obj);
options = check_tStart_tFinal(obj, options, checkName);

if ~hyb
    check_timeStep(options, obj, checkName);
    options = check_inputSet(obj, options);
end

% algorithm properties:

check_taylorTerms(options, obj);
check_zonotopeOrder(options, obj);
check_intermediateOrder(options, obj);



% OPTIONAL PARAMETERS -----------------------------------------------------

check_saveOrder(options, obj);



% WARNINGS ----------------------------------------------------------------

% do not print warnings if called from HA obj
if ~suppressWarnings

    % warnings for unused options (overdefined)
    validFields = getValidFields(class(obj));
    warning(printRedundancies(options,validFields));
    
end


% INTERNAL SETTINGS -------------------------------------------------------

if ~hyb
    options = set_inputSet(obj, options, checkName);
end

end

%------------- END OF CODE --------------
