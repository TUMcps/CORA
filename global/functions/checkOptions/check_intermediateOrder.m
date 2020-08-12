function check_intermediateOrder(options, obj)
% check_intermediateOrder - checks if options.intermediateOrder
%  1) exists
%  2) take an allowed values
%
% Syntax:
%    check_intermediateOrder(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - object of system for case differentiation
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

strct = 'options';
option = 'intermediateOrder';
if isa(obj, 'linParamSys')
    % linParamSys
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,strct));
    elseif isfield(options,option)
        if options.intermediateOrder < 1
            error(printOptionOutOfRange(obj,option,strct));
        end
    end
else
    % nonlinSys, nonlinDASys, nonlinParamSys
    % intermediateOrder has to be >= 1
    % if alg = 'poly' -> intermediateOrder has to exist
    if strcmp(options.alg,'poly')
        if ~isfield(options,option)
            error(printOptionMissing(obj,option,strct));
        end
    elseif isfield(options,option)
        if options.intermediateOrder < 1
            error(printOptionOutOfRange(obj,option,strct));
        end
    end
end

end

%------------- END OF CODE --------------
