function options = checkOptionsReach(obj,options,hyb)
% checkOptionsReach - checks if all necessary options
%    are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options,hyb)
%
% Inputs:
%    obj     - linearSys object (unused)
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
% Written:      28-Jan-2019
% Last update:  17-Aug-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';
if nargin ~= 3
    hyb = 0;
end
suppressWarnings = hyb;

% DEFAULT OPTIONS ---------------------------------------------------------

options = check_verbose(options, obj);
options = check_linAlg(options, obj);
options = check_reductionTechnique(options, obj);
options = check_tStart_tFinal(obj, options, checkName);
% only options.linAlg = 'adap'
options = check_error(options, obj);

% MANDATORY OPTIONS -------------------------------------------------------

% no timeStep, taylorTerms and zonotopeOrder in adaptive algorithm
if ~strcmp(options.linAlg,'adap')
    if ~hyb
        % note: tStart needed before check_timeStep called
        check_timeStep(options, obj, checkName);
    end

    % algorithm properties:
    check_taylorTerms(options, obj);
    check_zonotopeOrder(options, obj);
end

% input sizes must match
% skip if called from hybrid (uLoc etc. there)
if ~hyb
    check_R0(options, obj);
    options = check_inputSet(obj, options);
end


% OPTIONAL PARAMETERS -----------------------------------------------------

options = check_specification(obj, options);
check_saveOrder(options, obj);
% only options.linAlg = 'decomp':
check_partition(options, obj);
% only options.linAlg = 'krylov':
check_krylovError(options, obj);


% WARNINGS ----------------------------------------------------------------

% do not print warnings if called from HA obj
if ~suppressWarnings

    % warnings for unused options (overdefined)
    validFields = getValidFields('linearSys','all',options.linAlg);
    warning(printRedundancies(options,validFields));
    
end


% INTERNAL SETTINGS -------------------------------------------------------

if ~hyb
    options = set_inputSet(obj, options, checkName);
end


end

%------------- END OF CODE --------------
