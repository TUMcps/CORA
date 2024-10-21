function cZ = conZonotope(SpS)
% conZonotope - overapproximates a spectrahedron by a constrained zonotope
%
% Syntax:
%    cZ = conZonotope(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    cZ - zonotope object
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    cZ = conZonotope(SpS);
%
%    figure; hold on;
%    plot(SpS,[1,2],'b');
%    plot(cZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       04-August-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cZ = conZonotope(zonotope(SpS));

% ------------------------------ END OF CODE ------------------------------
