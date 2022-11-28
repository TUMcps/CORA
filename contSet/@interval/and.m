function res = and(I,S)
% and - computes intersection, overloades '&' operator of intervals
%
% Syntax:  
%    res = and(I,S)
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
% See also: zonotope/and

% Author:       Matthias Althoff
% Written:      26-June-2015
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_and('interval',I,S);

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % potential re-ordering
    I = vars{1}; S = vars{2};
end


% different cases depending on the class of the summand
if isa(S,'interval')

    % compute intersection
    lb = max(I.inf, S.inf);
    ub = min(I.sup, S.sup);

    % check if result is empty
    tmp = lb - ub;
    if all(all(tmp < eps | withinTol(tmp,eps)))
        res = interval(min([lb,ub],[],2),max([lb,ub],[],2));
    else
        res = [];
    end

elseif isa(S,'halfspace') || isa(S,'conHyperplane')

    % convert to conZonotope
    cZ = conZonotope(I);

    % compute intersection
    res = cZ & S;

    % ecnlose intersection by interval
    res = interval(res);

elseif isa(S,'levelSet')

    res = S & I;

elseif isa(S,'zonotope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'mptPolytope') || ...
       isa(S,'conPolyZono')

    res = S & I;

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',I,S));
    
end

    
end

%------------- END OF CODE --------------