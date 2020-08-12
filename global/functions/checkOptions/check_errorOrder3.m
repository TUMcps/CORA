function check_errorOrder3(options, obj)
% check_errorOrder3 - checks if options.errorOrder3
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_errorOrder3(options, obj)
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
% See also: check_errorOrder
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      21-April-2020
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'errorOrder3';
% errorOrder3 has to be >= 1
if strcmp(options.alg,'poly') && options.tensorOrder >= 4
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,'options'));
    elseif options.errorOrder3 < 1
        error(printOptionOutOfRange(obj,option));
    end
end

end

%------------- END OF CODE --------------
