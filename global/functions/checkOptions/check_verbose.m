function options = check_verbose(options, obj)
% check_verbose - checks if options.verbose
%  1) takes an allowed value
%
% Syntax:
%    options = check_verbose(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    options - options for object
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
% Last update:  26-Aug-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

% verbose has to be 0 or 1
option = 'verbose';
strct = 'options';
defValue = getDefaultOption(option);
% check verbose ... optional, default at 0
if ~isfield(options, option)
    options.verbose = false;
elseif options.verbose == defValue
    % not necessary since default value
    printDefaultValue(obj,option,defValue);
elseif ~(options.verbose == 0 || options.verbose == 1)
    error(printOptionOutOfRange(obj, option, strct));
end


end

%------------- END OF CODE --------------

