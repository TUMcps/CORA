function res = and_(I,S,varargin)
% and_ - computes intersection, overloades '&' operator of STL intervals
%
% Syntax:
%    res = and_(I,S)
%
% Inputs:
%    I - stlInterval object
%    S - contSet object
%
% Outputs:
%    res - intersection of interval objects
%
% Example: 
%    I1 = stlInterval(0,3);
%    I2 = stlInterval(1,2,false,true);
%    res = I1 & I2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Intersection operator for intervals
if isa(S,'stlInterval')
    if isemptyobject(I) || isemptyobject(S)
        res = stlInterval();
        return;
    end

    lb = max(I.lower,S.lower);
    ub = min(I.upper,S.upper);
    if I.lower > S.lower
        lc = I.leftClosed;
    elseif I.lower < S.lower
        lc = S.leftClosed;
    else
        lc = I.leftClosed && S.leftClosed;
    end
    if I.upper < S.upper
        rc = I.rightClosed;
    elseif I.upper > S.upper
        rc = S.rightClosed;
    else
        rc = I.rightClosed && S.rightClosed;
    end
    
    res = stlInterval(lb,ub,lc,rc);
else
    throw(CORAerror('CORA:noops',I,S));
end

% ------------------------------ END OF CODE ------------------------------
