function SpS = spectraShadow(C)
% spectraShadow - Converts a capsule to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    S - spectraShadow object
%
% Example:
%    c = [1;2];
%    g = [2;1];
%    r = 1;
%    C = capsule(c,g,r);
%    S = spectraShadow(C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       03-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = C.c; % center
g = C.g; % generator
r = C.r; % radius

n = dim(C);

% Create the stretched sphere
sphere = r * ellipsoid(eye(n), zeros([n 1]));
sphere = spectraShadow(sphere);

% Now create the segment, as a zonotope
segment = zonotope(c, g);
segment = spectraShadow(segment);

% Spectrahedra know how to deal with Minkowski additions, so let them do
% the heavy lifting

SpS = sphere + segment;

% Additional properties
SpS.bounded.val = true;
SpS.emptySet.val = representsa_(C,'emptySet',1e-10);
SpS.fullDim.val = isFullDim(C);
SpS.center.val = center(C);

% ------------------------------ END OF CODE ------------------------------
