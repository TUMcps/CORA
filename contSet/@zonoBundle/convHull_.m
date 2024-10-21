function S_out = convHull_(zB,S,varargin)
% convHull_ - computes the convex hull of two zonotope bundles
%
% Syntax:
%    res = convHull_(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object or numeric
%
% Outputs:
%    S_out - zonoBundle enclosing the convex hull
%
% Example: 
%    I = interval([4;3],[6;6]);
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
% 
%    zB_ = convHull(zB,I);
% 
%    figure; hold on;
%    plot(I);
%    plot(zB);
%    plot(zB_);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull, conZonotope/convHull_

% Authors:       Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% zonotope bundles are convex
if nargin == 1
    S_out = zB;
    return;
end

% ensure that numeric is second input argument
[zB,S] = reorderNumeric(zB,S);

% check dimensions
equalDimCheck(zB,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < zB.precedence
    S_out = convHull(S,zB,varargin{:});
    return
end

% check for empty case
if representsa_(zB,'emptySet',1e-12)
    S_out = S;
    return
elseif representsa_(S,'emptySet',1e-12)
    S_out = zB;
    return
end

% general idea: convert to polytope and back
try
    P2 = polytope(S);
catch ME
    throw(CORAerror('CORA:noops',zB,S));
end
P1 = polytope(zB);

% compute convex hull
P = convHull(P1,P2);

% convert back to a zonotope bundle
S_out = zonoBundle(P);

% ------------------------------ END OF CODE ------------------------------
