function check_x0(options, obj)
% check_x0 - checks if options.verbose
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
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      03-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

option = 'x0';
% size of x0 must match
if ~isfield(options,option)
    error(printOptionMissing(obj,option,'params'));
elseif length(options.x0) ~= size(center(options.R0),1)
    error(printOptionSpecificError(obj,option,...
            'R0 and x0 must have the same dimension.'));
end


end

%------------- END OF CODE --------------

