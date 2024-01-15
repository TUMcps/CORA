function res = contains_(I,S,type,tol,varargin)
% contains_ - determines if an interval contains a set or a point
%
% Syntax:
%    res = contains_(I,S)
%    res = contains_(I,S,type)
%    res = contains_(I,S,type,tol)
%
% Inputs:
%    I - interval object
%    S - contSet object or single point
%    type - 'exact' or 'approx'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([-1;-2],[2;3]);
%    I2 = interval([0;0],[1;1]);
%
%    contains(I1,I2)
%    contains(I1,[1;2])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper
% Written:       16-May-2018
% Last update:   02-September-2019
%                15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% init result
res = false;

% set in empty set
if representsa_(I,'emptySet',0)
    res = representsa_(S,'emptySet',0);
    return
end

% point in interval containment
if isnumeric(S)
    
    res = all( (I.inf < S | withinTol(I.inf,S,tol)) ...
            & (I.sup > S | withinTol(I.sup,S,tol)) , 1);

% interval in interval containment
elseif isa(S,'interval')

    % TODO: use withinTol(.,.,tol) below
    if all(I.sup >= S.sup) && all(I.inf <= S.inf)
        res = true;
    end

% other set in interval containment
else

    res = contains(polytope(I),S,type,tol);

end

% ------------------------------ END OF CODE ------------------------------
