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
%                28-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence (all but interval-interval case)
if isa(S,'contSet') && S.precedence < I.precedence
    res = and_(S,I,varargin{:});
    return
end

% interval-interval case

% compute intersection
lb = max(I.inf, S.inf);
ub = min(I.sup, S.sup);

% check if intersection is empty
tmp = lb - ub;
if all(tmp <= eps, 'all')
    res = interval(min([lb,ub],[],2),max([lb,ub],[],2));
else
    res = interval.empty(prod(dim(S)));
end


% elseif isa(S,'polytope') && (representsa_(S,'halfspace',1e-12) || representsa_(S,'conHyperplane',1e-12))
% 
%     % convert to conZonotope
%     cZ = conZonotope(I);
% 
%     % compute intersection
%     res = and_(cZ,S,'exact');
% 
%     % ecnlose intersection by interval
%     res = interval(res);
% 
% elseif isa(S,'zonotope') || isa(S,'conZonotope') || ...
%        isa(S,'zonoBundle') || isa(S,'polytope')
% 
%     res = and_(S,I,'exact');
% 
% end

% ------------------------------ END OF CODE ------------------------------
