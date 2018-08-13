function res = prod( obj, dim )
% prod - product of array elements
%
% Syntax:  
%    res = prod( obj, dim )
%
% Inputs:
%    obj - input array (interval)
%    dim - 1 - product of column's elements; 2 - of row's elements
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

% Author:       Dmitry Grebenyuk
% Written:      24-October-2017
% Last update:  
% Last revision:---

    if dim == 1 % reduce to a row 
        S.type='()'; % to avoid Matlab's bug
        S.subs={1,':'};
        res = subsref(obj,S);
        for i = 2:size(obj, dim)
            S.subs={i,':'};
            res = res .* subsref(obj,S);
        end
    elseif dim == 2 % reduce to a column
        S.type='()';
        S.subs={':', 1};
        res = subsref(obj,S);
        for i = 2:size(obj, dim)
            S.subs={':', i};
            res = res .* subsref(obj,S);
        end
    else
        error ('Wrong input')
    end
end

