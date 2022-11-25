function res = cartProd(Int1,Int2)
% cartProd - returns the Cartesian product of two intervals
%
% Syntax:  
%    res = cartProd(Int1,Int2)
%
% Inputs:
%    Int1 - interval object
%    Int2 - interval object
%
% Outputs:
%    res - Cartesian product of Int1 and Int2
%
% Example: 
%    a = interval([-1;-3],[1;6]);
%    b = interval([-2;1], [4;2]);
%    c = cartProd(a,b)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartProd

% Author:       Mark Wetzlinger
% Written:      18-Sep-2019
% Last update:  24-Sep-2019
%               05-May-2020 (standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% first or second set is interval
if isa(Int1,'interval')

    % different cases for different set representations
    if isa(Int2,'interval')

        if isempty(Int1)
            res = Int2;
        elseif isempty(Int2)
            res = Int1;
        elseif size(Int1,2) == 1 && size(Int2,2) == 1 % vertcat
            res = interval([Int1.inf;Int2.inf],[Int1.sup;Int2.sup]);
        elseif size(Int1,1) == 1 && size(Int2,1) == 1 % horzcat
            res = interval([Int1.inf,Int2.inf],[Int1.sup,Int2.sup]);
        else
            error("No Cartesian product possible");
        end

    elseif isnumeric(Int2)

        res = interval([Int1.inf;Int2],[Int1.sup;Int2]);

    % different cases for different set representations
    elseif isa(Int2,'zonotope') 
        res = cartProd(zonotope(Int1),Int2);
    elseif isa(Int2,'conZonotope')
        res = cartProd(conZonotope(Int1),Int2);
    elseif isa(Int2,'zonoBundle')
        res = cartProd(zonoBundle(Int1),Int2);
    elseif isa(Int2,'mptPolytope')
        res = cartProd(mptPolytope(Int1),Int2);
    elseif isa(Int2,'polyZonotope')
        res = cartProd(polyZonotope(Int1),Int2);
    elseif isa(Int2,'conPolyZono')
        res = cartProd(conPolyZono(Int1),Int2);
    else
        % throw error for given arguments
        error(noops(Int1,Int2));
    end

else

    if isnumeric(Int1)
        res = interval([Int1;Int2.inf],[Int1;Int2.sup]);
    else
        % throw error for given arguments
        error(noops(Int1,Int2));
    end  
end

    
end

%------------- END OF CODE --------------