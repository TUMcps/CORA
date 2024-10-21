function S_out = convHull_(Z,S,varargin)
% convHull_ - computes an enclosure for the convex hull of a zonotope and
%    another set or a point
%
% Syntax:
%    Z = convHull_(Z,S)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%
% Outputs:
%    Z - zonotope enclosing the convex hull
%
% Example: 
%    Z1 = zonotope([2 1 0; 2 0 1]);
%    Z2 = zonotope([-2 1 0; -2 0 1]);
%
%    Z = convHull(Z1,Z2);
%
%    figure; hold on;
%    plot(Z1,[1,2],'FaceColor','r');
%    plot(Z2,[1,2],'FaceColor','b');
%    plot(Z,[1,2],'g','LineWidth',3);
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

% zonotope is already convex
if nargin == 1
    S_out = Z;
    return;
end

% ensure that numeric is second input argument
[Z,S] = reorderNumeric(Z,S);

% check dimensions
equalDimCheck(Z,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < Z.precedence
    S_out = convHull(S,Z,varargin{:});
    return
end

% convex hull with empty set
if representsa_(S,'emptySet',eps)
    S_out = Z;
    return
elseif representsa_(Z,'emptySet',eps)
    S_out = S;
    return
end

% use enclose method
S_out = enclose(Z,zonotope(S));

% ------------------------------ END OF CODE ------------------------------
