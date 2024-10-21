function S_out = convHull_(C,S,varargin)
% convHull_ - computes an enclosure for the convex hull of a zonotope and
%    another set or a point
%
% Syntax:
%    C = convHull_(C)
%    C = convHull_(C,S)
%
% Inputs:
%    C - capsule object
%    S - contSet object, numerical vector
%
% Outputs:
%    S_out - convex hull
%
% Example: 
%    C1 = capsule([1;0],[1;1],0.5);
%    C2 = capsule([-2;-1],[-1;1],0.2);
%    C = convHull(C1,C2);
% 
%    figure; hold on;
%    plot(C1,[1,2],'FaceColor','r');
%    plot(C2,[1,2],'FaceColor','b');
%    plot(C,[1,2],'g','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull

% Authors:       Mark Wetzlinger
% Written:       17-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% capsule is already convex
if nargin == 1
    S_out = C;
    return;
end

% ensure that numeric is second input argument
[C,S] = reorderNumeric(C,S);

% check dimensions
equalDimCheck(C,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < C.precedence
    S_out = convHull(S,C,varargin{:});
    return
end

% convex hull with empty set is the other set
if representsa_(S,'emptySet',eps)
    return;
end

% capsule-capsule case
if isa(S,'capsule')

    % choose new center as mid-point between two centers
    c_ = (C.c + S.c)/2;

    % generator as distance from new center to one of the old centers
    g_ = C.c - c_;

    % radius to enclose both capsules
    r_ = max([vecnorm(C.g) + C.r,vecnorm(S.g) + S.r]);

    % instantiate result
    S_out = capsule(c_,g_,r_);
    return
end

% throw error for given arguments
throw(CORAerror('CORA:noops',C,S));

% ------------------------------ END OF CODE ------------------------------
