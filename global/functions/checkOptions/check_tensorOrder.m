function check_tensorOrder(options, obj)
% check_tensorOrder - checks if options.tensorOrder
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_tensorOrder(options, obj)
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
% Last update:  21-April-2020 (split in lin/poly)
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'tensorOrder';
strct = 'options';
% tensorOrder has to be integer and
%  at least 2 if options.alg = lin
%  at least 3 if options.alg = poly
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
else
    maxOrder = inf;
    if isa(obj,'nonlinearSysDT')
        minOrder = 2;
    elseif strcmp(options.alg,'lin')
        minOrder = 2;
    elseif strcmp(options.alg,'linRem')
        minOrder = 2;
        maxOrder = 2;
    elseif strcmp(options.alg,'poly')
        minOrder = 3;
    end
    if options.tensorOrder < minOrder || options.tensorOrder > maxOrder || mod(options.tensorOrder,1.0) ~= 0
        error(printOptionOutOfRange(obj,option,strct));
    end
end

end

%------------- END OF CODE --------------
