function options = checkOptionsReach(obj,options,hyb)
% checkOptionsReach - checks if all necessary options
%   are there and have valid values
%
% Syntax:  
%    options = checkOptionsReach(obj,options,hyb)
%
% Inputs:
%    obj     - linearSysDT object (unused)
%    options - options for linearSysDT (1x1 struct)
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      20-Mar-2020
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
options = check_reductionTechnique(options, obj);
options = check_tStart_tFinal(obj, options, checkName);

t = options.tStart:obj.dt:options.tFinal;
if t(end) ~= options.tFinal
   error('Final time has to be a multiple of the sampling time!'); 
end

% MANDATORY OPTIONS -------------------------------------------------------

% algorithm properties:
check_zonotopeOrder(options, obj);


% input sizes must match
% skip if called from hybrid (uLoc etc. there)
if ~hyb
    check_R0(options, obj);
    options = check_inputSet(obj, options);
end


% OPTIONAL PARAMETERS -----------------------------------------------------

options = check_specification(obj, options);
check_saveOrder(options, obj);


% WARNINGS ----------------------------------------------------------------

% do not print warnings if called from HA obj
if ~suppressWarnings

    % warnings for unused options (overdefined)
    validFields = getValidFields('linearSysDT','all');
    warning(printRedundancies(options,validFields));
    
end


% INTERNAL SETTINGS -------------------------------------------------------

if ~hyb
    options = set_inputSet(obj, options, checkName);
end

end

%------------- END OF CODE --------------
