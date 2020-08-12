function options = checkOptionsReach(obj,options,hyb)
% checkOptionsReach - checks if all necessary options
%   are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options,hyb)
%
% Inputs:
%    obj     - linearSys object (unused)
%    options - options for linearSys (1x1 struct)
%    hyb     - if called from hybrid class
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
% Written:      05-Feb-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsReach';
if nargin == 3
    suppressWarnings = hyb;
else
    suppressWarnings = 0;
end

% DEFAULT OPTIONS ---------------------------------------------------------

options = check_tStart_tFinal(obj, options, checkName);
options = check_reductionTechnique(options, obj);
options = check_verbose(options, obj);
check_lagrangeRem(options, obj);

t = options.tStart:obj.dt:options.tFinal;
if t(end) ~= options.tFinal
   error('Final time has to be a multiple of the sampling time!'); 
end

% MANDATORY OPTIONS -------------------------------------------------------

check_R0(options, obj);
options = check_inputSet(obj, options);

% algorithm properties:

check_zonotopeOrder(options, obj);

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

