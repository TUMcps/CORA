function options = check_linAlg(options, obj)
% check_linAlg - checks if options.linAlg
%  1) takes an allowed value
%
% Syntax:
%    check_linAlg(options, obj)
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
% Written:      20-Dec-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'linAlg';
strct = 'options';

defValue = getDefaultOption(option);
% the chosen algorithm has to match a predefined string
validAlg = {'standard';'wrapping-free';'fromStart';'decomp';'krylov';'adap'};
if ~isfield(options,option)
    options.linAlg = defValue;
elseif strcmp(options.linAlg,defValue)
    % not necessary since default value
    printDefaultValue(obj,option,defValue);
elseif ~any(strcmp(validAlg,options.linAlg))
    error(printOptionOutOfRange(obj,option,strct));
end

end

%------------- END OF CODE --------------

