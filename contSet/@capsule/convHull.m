function C = convHull(C,S)
% convHull - computes an enclosure for the convex hull of a zonotope and
%    another set or a point
%
% Syntax:
%    C = convHull(C,S)
%
% Inputs:
%    C - capsule object
%    S - contSet object, numerical vector
%
% Outputs:
%    C - capsule enclosing the convex hull
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
% See also: conZonotope/convHull

% Authors:       Mark Wetzlinger
% Written:       17-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 1
    return;
end 
if representsa_(S,'emptySet',eps)
    return;
end

% determine zonotope object
[C,S] = findClassArg(C,S,'capsule');

% different cases depending on the class of the second summand
if isa(S,'capsule')

    % choose new center as mid-point between two centers
    c_ = (C.c + S.c)/2;

    % generator as distance from new center to one of the old centers
    g_ = C.c - c_;

    % radius to enclose both capsules
    r_ = max([vecnorm(C.g) + C.r,vecnorm(S.g) + S.r]);

    % instantiate result
    C = capsule(c_,g_,r_);

else

    % throw error for given arguments
    throw(CORAerror('CORA:noops',C,S));

end

% ------------------------------ END OF CODE ------------------------------
