function cPZ = conPolyZono(SpS)
% conPolyZono - overapproximates a spectrahedron by a constrained 
%    polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    A0 = eye(4);
%    A1 = blkdiag([1 0;0 -1],zeros(2));
%    A2 = blkdiag(zeros(2),[1 0;0 -1]);
%    SpS = spectraShadow([A0 A1 A2]);
%    cPZ = conPolyZono(SpS);
%
%    figure; hold on;
%    plot(SpS,[1,2],'b');
%    plot(cPZ,[1,2],'r');
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

cPZ = conPolyZono(zonotope(SpS));

% ------------------------------ END OF CODE ------------------------------
