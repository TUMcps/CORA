function options = checkOptionsSimulate(obj,options,cDsimRand)
% checkOptionsSimulate - checks if all necessary options
%    are there and have valid values
%
% Syntax:  
%    options = checkOptionsSimulate(obj,options,cDsimRand)
%
% Inputs:
%    obj       - linearSys object (unused)
%    options   - options for linearSys (1x1 struct)
%    cDsimRand - (optional) if called from @contDynamics > simulateRandom
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
% Written:      02-Feb-2019
% Last update:  18-Dec-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

checkName = 'checkOptionsSimulate';
if nargin < 3
    cDsimRand = 0;
end


% MANDATORY OPTIONS -------------------------------------------------------

check_metasim(options, obj, cDsimRand);

% other options
options = check_tStart_tFinal(obj, options, checkName);
check_timeStep(options, obj, checkName);
check_R0(options, obj);
options = check_inputSet(obj, options);
% only for nonlinDASys objects
check_y0guess(options, obj);


% WARNINGS ----------------------------------------------------------------

% warning, if both uTrans and uTransVec defined
if isfield(options,'uTrans') && isfield(options,'uTransVec')
    warning('Only uTransVec used, uTrans is therefore redundant');
end

% necessary fields
if cDsimRand
    method = 'simRand';
else
    method = 'RRT';
end
validFields = getValidFields('contDyn_sim','all',method);

% add y0guess if nonlinDASys object
if isa(obj, 'nonlinDASys')
    validFields = [validFields; 'y0guess'];
end

warning(printRedundancies(options,validFields));


% INTERNAL SETTINGS -------------------------------------------------------

options = set_inputSet(obj, options, checkName);

end

%------------- END OF CODE --------------
