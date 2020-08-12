function check_simplify(options, obj)
% check_simplify - checks if options.simplify
%  1) takes an allowed value
%
% Syntax:
%    check_simplify(options, obj)
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
% Last update:  14-Aug-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'options.lagrangeRem';
option = 'simplify';
% simplify has to be one of the strings below
validStrings = {'none';'simplify';'collect';'optimize'};
if isfield(options,option)
    if ~any(strcmp(validStrings,options.simplify))
        error(printOptionOutOfRange(obj,option,strct));
    end
end

end

%------------- END OF CODE --------------

