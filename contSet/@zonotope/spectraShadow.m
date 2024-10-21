function SpS = spectraShadow(Z)
% spectraShadow - Converts a zonotope to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    SpS - spectraShadow object
%
% Example:
%    Z = zonotope([0;0],eye(2));
%    SpS = spectraShadow(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% The conZonotope implementation is more general, so do it that way
SpS = spectraShadow(conZonotope(Z));

% Additional properties
SpS.bounded.val = true;
SpS.emptySet.val = representsa_(Z,'emptySet',1e-10);
SpS.fullDim.val = isFullDim(Z);
SpS.center.val = center(Z);

% ------------------------------ END OF CODE ------------------------------
