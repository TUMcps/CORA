function check_pha_x0(options, obj)
% check_pha_x0 - checks if options.x0
%  1) takes an allowed value
%
% Syntax:
%    check_x0(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      14-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    option = 'x0';
    
    % size of x0 must match
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,'params'));
    elseif ~all(size(options.x0) == [obj.numStates,1])
        error(printOptionOutOfRange(obj,option,'params.'));
    end
end

%------------- END OF CODE --------------