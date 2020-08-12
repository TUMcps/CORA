function options  = check_reductionInterval(options, obj)
% check_reductionInterval - checks if options.reductionInterval
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_reductionInterval(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - checkOptions* function (for error tracing)
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

strct = 'options';
option = 'reductionInterval';
defValue = getDefaultOption(option);
% reductionInterval has to be integer and at least 0
if ~isfield(options,option)
    options.reductionInterval = defValue;
elseif options.reductionTechnique == defValue
    % not necessary since default value
    printDefaultValue(obj,option,defValue);
elseif options.reductionInterval < 0 || mod(options.reductionInterval,1.0) ~= 0
    if ~isinf(options.reductionInterval)
        error(printOptionOutOfRange(obj,option,strct));
    end
end

end

%------------- END OF CODE --------------

