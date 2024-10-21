function res = isIntersecting_(Z,S,type,tol,varargin)
% isIntersecting_ - determines if zonotope intersects a set
%
% Syntax:
%    res = isIntersecting_(Z,S,type,tol)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope([0 1 1 0;0 1 0 1]);
%    Z2 = zonotope([2 -1 1 0;2 1 0 1]);
%    Z3 = zonotope([3.5 -1 1 0;3 1 0 1]);
% 
%    isIntersecting(Z1,Z2)
%    isIntersecting(Z1,Z3)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, conZonotope/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       21-November-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[Z,S] = reorderNumeric(Z,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < Z.precedence
    res = isIntersecting_(S,Z,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(Z,S,type,tol);
    return
end

% sets must not be empty
if representsa_(Z,'emptySet',0) || representsa_(S,'emptySet',0)
    res = false;
    return
end

% general idea: convert to constrained zonotope and check for intersection
% unless S is unbounded as unbounded sets are not convertible to conZonotopes
if isa(S,'contSet')
    if isBounded(S) 
        % both bounded, convert to constrained zonotope
        res = isIntersecting_(conZonotope(Z),conZonotope(S),type,tol,varargin{:});
    else
        % S unbounded, use polytope function
        res = isIntersecting_(polytope(S),Z,type,tol,varargin{:});
    end
    return
end

throw(CORAerror('CORA:noops',Z,S));

% ------------------------------ END OF CODE ------------------------------
