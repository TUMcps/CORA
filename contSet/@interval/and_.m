function res = and_(I,S,varargin)
% and_ - computes intersection, overloades '&' operator of intervals
%
% Syntax:
%    res = and_(I,S)
%
% Inputs:
%    I - interval object
%    S - contSet object
%
% Outputs:
%    res - intersection of interval objects
%
% Example: 
%    I1 = interval([1; -1], [2; 1]);
%    I2 = interval([1.5; -2], [2.5; 0]);
%    res = I1 & I2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, zonotope/and_

% Authors:       Matthias Althoff
% Written:       26-June-2015
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

% different cases depending on the class of the summand
if isa(S,'interval')

    % compute intersection
    lb = max(I.inf, S.inf);
    ub = min(I.sup, S.sup);

    % check if result is empty
    tmp = lb - ub;
    if all(tmp <= eps, 'all')
        res = interval(min([lb,ub],[],2),max([lb,ub],[],2));
    else
        res = interval.empty(dim(S));
    end

elseif isa(S,'halfspace') || isa(S,'conHyperplane')

    % convert to conZonotope
    cZ = conZonotope(I);

    % compute intersection
    res = and_(cZ,S,'exact');

    % ecnlose intersection by interval
    res = interval(res);

elseif isa(S,'levelSet')

    res = and_(S,I,'exact');

elseif isa(S,'zonotope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'polytope') || ...
       isa(S,'conPolyZono')

    res = and_(S,I,'exact');

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',I,S));
    
end

% ------------------------------ END OF CODE ------------------------------
