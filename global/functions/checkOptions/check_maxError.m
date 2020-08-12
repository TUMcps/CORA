function options = check_maxError(options, obj)
% check_maxError - checks if options.maxError
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    options = check_maxError(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - updated options for object
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
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'maxError';
defValue = getDefaultOption(option);
% maxError same dim as R0, every entry >= 0
if ~isfield(options,option)
    dim = options.R0.dim;
    options.maxError = defValue * ones(dim,1);
elseif size(center(options.R0),1) ~= length(options.maxError)
    error(printOptionSpecificError(obj,option,...
            'R0 and maxError must have the same dimension.'));
elseif any((options.maxError >= 0) == 0)
    error(printOptionSpecificError(obj,option,...
            'Every entry in maxError has to be >= 0.'));
end

end

%------------- END OF CODE --------------

