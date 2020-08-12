function res = prod(obj,n)
% prod - product of array elements
%
% Syntax:  
%    res = prod(obj,n)
%
% Inputs:
%    obj - interval object
%    n - dimension:
%           1 - product of column's elements
%           2 - product of row's elements
%
% Outputs:
%    res - interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Dmitry Grebenyuk
% Written:      24-October-2017
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------

if n == 1 % reduce to a row 
    S.type='()'; % to avoid Matlab's bug
    S.subs={1,':'};
    res = subsref(obj,S);
    for i = 2:size(obj, n)
        S.subs={i,':'};
        res = res .* subsref(obj,S);
    end
elseif n == 2 % reduce to a column
    S.type='()';
    S.subs={':', 1};
    res = subsref(obj,S);
    for i = 2:size(obj, n)
        S.subs={':', i};
        res = res .* subsref(obj,S);
    end
else
    error ('Wrong input')
end

%------------- END OF CODE --------------
