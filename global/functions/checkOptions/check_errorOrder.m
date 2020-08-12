function check_errorOrder(options, obj)
% check_errorOrder - checks if options.errorOrder
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_errorOrder(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - object of system for case differentiation
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
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'errorOrder';
strct = 'options';
if isa(obj, 'nonlinearSysDT')
    % nonlinearSysDT
    % errorOrder has to be >= 1
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif options.errorOrder < 1
        error(printOptionOutOfRange(obj,option,strct));
    end
else
    % nonlinearSys, nonlinDASys, nonlinParamSys
    % errorOrder has to be >= 1
    % if alg = 'poly' -> errorOrder has to exist
    if strcmp(options.alg,'poly') || ...
        (strcmp(options.alg,'lin') && options.tensorOrder > 2)
        if ~isfield(options,option)
            error(printOptionMissing(obj,option,strct));
        end
    elseif isfield(options,option)
        if options.errorOrder < 1
            error(printOptionOutOfRange(obj,option,strct));
        end
    end
end

end

%------------- END OF CODE --------------
