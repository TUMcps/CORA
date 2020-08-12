function options = checkOptionsReach(obj,options,hyb)
% checkOptionsReach - checks if all necessary options
%    are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options,hyb)
%
% Inputs:
%    obj     - nonlinDASys object (unused)
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
% Written:      17-Feb-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';
if nargin == 3
    suppressWarnings = hyb;
else
    hyb = 0;
    suppressWarnings = hyb;
end

% DEFAULT OPTIONS ---------------------------------------------------------

options = check_verbose(options, obj);
options = check_reductionTechnique(options, obj);
options = check_linAlg(options, obj);
options = check_reductionInterval(options, obj);

if isfield(options,'alg')
    warning('The following params/options have been set, but are redundant: alg')
end
options.alg = 'lin'; % for checks only

% MANDATORY OPTIONS -------------------------------------------------------

check_R0(options, obj);
check_y0guess(options, obj);

options = check_tStart_tFinal(obj, options, checkName);
check_timeStep(options, obj, checkName);

options = check_inputSet(obj, options);

% algorithm properties:
check_taylorTerms(options, obj);
check_zonotopeOrder(options, obj);
check_tensorOrder(options, obj);
check_maxError_xy(options, obj);


% OPTIONAL PARAMETERS -----------------------------------------------------

options = check_maxError(options, obj);
check_intermediateOrder(options, obj);
check_errorOrder(options, obj);
check_errorOrder3(options, obj);
check_replacements(options, obj);
check_tensorParallel(options, obj);
check_lagrangeRem(options, obj);
check_simplify(options, obj);
check_saveOrder(options, obj);


% WARNINGS ----------------------------------------------------------------

options = rmfield(options,'alg');

% do not print warnings if called from HA obj
if ~suppressWarnings

    % warnings for unused options (overdefined)
    validFields = getValidFields(class(obj));
    warning(printRedundancies(options,validFields));
    
end


% INTERNAL SETTINGS -------------------------------------------------------

if ~suppressWarnings
    options = set_inputSet(obj, options, checkName);
end


end

%------------- END OF CODE --------------

