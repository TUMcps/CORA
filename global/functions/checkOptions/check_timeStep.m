function check_timeStep(options, obj, checkName)
% check_timeStep - checks if options.timeStep
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_timeStep(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - system object
%    checkName - if reach or simulate checked
%
% Outputs:
%    -
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
% Written:      04-Mar-2019
% Last update:  17-Aug-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'timeStep';
strct = 'options';
% timeStep is redundant if called from @contDynamics > simulateRandom,
%   but mandatory in reach (except linearSys 'adap')
if strcmp(checkName,'checkOptionsReach') && ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
elseif isfield(options,option) && ...
    (options.timeStep > (options.tFinal - options.tStart) + 10*eps || ...
        options.timeStep <= 0)
    error(printOptionOutOfRange(obj,option,strct));
end



end

%------------- END OF CODE --------------
