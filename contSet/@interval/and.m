function res = and(Int1,Int2)
% and - computes intersection of intervals.
%       Overloades '&' operator for intervals.
%
% Syntax:  
%    res = and(Int1,Int2)
%
% Inputs:
%    Int1 - first interval object
%    Int2 - second interval object
%
% Outputs:
%    res - resulting interval object
%
% Example: 
%    a = interval([1;-1], [2; 1]);
%    b = interval([1.5; -2], [2.5; 0]);
%    c = a & b
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/and

% Author:       Matthias Althoff
% Written:      26-June-2015
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% determine the interval object
if ~isa(Int1,'interval')
    temp = Int1;
    Int1 = Int2;
    Int2 = temp;
end

% different cases depending on the class of the summand
if isa(Int2,'interval')

    % compute intersection
    inf = max(Int1.inf, Int2.inf);
    sup = min(Int1.sup, Int2.sup);

    % check if result is empty
    if all(all(inf - sup <= eps))
        res = interval(min([inf,sup],[],2),max([inf,sup],[],2));
    else
        res = [];
    end

elseif isa(Int2,'halfspace') || isa(Int2,'conHyperplane')

    % convert to conZonotope
    cZ = conZonotope(Int1);

    % compute intersection
    res = cZ & Int2;

    % ecnlose intersection by interval
    res = interval(res);

elseif isa(Int2,'levelSet')

    res = Int2 & Int1;

elseif isa(Int2,'zonotope') || isa(Int2,'conZonotope') || ...
       isa(Int2,'zonoBundle') || isa(Int2,'mptPolytope') || ...
       isa(Int2,'conPolyZono')

    res = Int2 & Int1;

else
    
    % throw error for given arguments
    error(noops(Int1,Int2));
    
end

    
end

%------------- END OF CODE --------------