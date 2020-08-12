function check_R0(options, obj)
% check_R0 - checks if options.R0
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_R0(options, obj)
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
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      04-Mar-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'params';
option = 'R0';
if isa(obj, 'parallelHybridAutomaton')
    % size of R0 must be equal to total number of states
    totalStates = 0;
    for i=1:length(obj.bindsStates)
        totalStates = totalStates + length(obj.bindsStates{i});
    end
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif size(center(options.R0),1) ~= totalStates
        error(printOptionOutOfRange(obj,option,strct));
    end
else
    % size of R0 must have same dimensionality as system
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif size(center(options.R0),1) ~= obj.dim
        error(printOptionSpecificError(obj,option,...
            'R0 and object must have the same dimension.'));
    end
end

end

%------------- END OF CODE --------------
