function res = isIntersecting_(zB,S,type,tol,varargin)
% isIntersecting_ - determines if zonotope bundle intersects a set
%
% Syntax:
%    res = isIntersecting_(zB,S,type,tol)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([2;2],[4;4]);
%    I2 = interval([3.5;3],[5;5]);
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%
%    isIntersecting(zB,I1)
%    isIntersecting(zB,I2)
%
%    figure; hold on
%    plot(zB,[1,2],'b');
%    plot(I1,[1,2],'g');
%    
%    figure; hold on
%    plot(zB,[1,2],'b');
%    plot(I2,[1,2],'r');
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
[zB,S] = reorderNumeric(zB,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < zB.precedence
    res = isIntersecting_(S,zB,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(zB,S,type,tol);
    return
end

% sets must not be empty (LP too costly...)
% if representsa_(zB,'emptySet',0) || representsa_(S,'emptySet',0)
%     res = false;
%     return
% end

if isa(S,'interval')
    res = isIntersecting_(polytope(S),zB,type,tol,varargin{:});
    return
end

% general idea: convert to constrained zonotope and check for intersection
if isa(S,'contSet')
    res = isIntersecting_(conZonotope(zB),conZonotope(S),type,tol,varargin{:});
    return
end

throw(CORAerror('CORA:noops',zB,S));

% ------------------------------ END OF CODE ------------------------------
