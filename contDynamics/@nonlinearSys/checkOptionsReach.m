function options = checkOptionsReach(obj,options,hyb)
% checkOptionsReach - checks if all necessary options
%    are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options,hyb)
%
% Inputs:
%    obj     - nonlinearSys object (unused)
%    options - options for nonlinearSys (1x1 struct)
%    hyb     - called from hybrid Automaton (0/1)
%
% Outputs:
%    options - options for nonlinearSys (1x1 struct)
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
% Written:      05-Feb-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';
if nargin ~= 3
    hyb = 0;
end
suppressWarnings = hyb;

% DEFAULT OPTIONS ---------------------------------------------------------

options = check_verbose(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
options = check_linAlg(options, obj);
options = check_reductionTechnique(options, obj);
options = check_maxError(options, obj);
options = check_reductionInterval(options, obj);
% using polyZonotope class...
options = check_polyZono(options, obj);


% MANDATORY OPTIONS -------------------------------------------------------

check_R0(options, obj);

if ~hyb
    % skip if called from hybrid (time step and input depending on loc)
    options = check_inputSet(obj, options);
    check_timeStep(options, obj, checkName);
end

% algorithm properties:

check_taylorTerms(options, obj);
check_zonotopeOrder(options, obj);
check_alg(options, obj);
check_tensorOrder(options, obj);

% OPTIONAL PARAMETERS -----------------------------------------------------

check_intermediateOrder(options, obj);
check_errorOrder(options, obj);
check_errorOrder3(options, obj);
check_lagrangeRem(options, obj);
check_saveOrder(options, obj);


% WARNINGS ----------------------------------------------------------------

% do not print warnings if called from HA obj
if ~suppressWarnings

    % warnings for unused options (overdefined)
    validFields = getValidFields(class(obj));
    warning(printRedundancies(options,validFields));
    
end


% INTERNAL SETTINGS -------------------------------------------------------

options = set_volApproxMethod(options);
if ~hyb
    options = set_inputSet(obj, options, checkName);
end

end

%------------- END OF CODE --------------

