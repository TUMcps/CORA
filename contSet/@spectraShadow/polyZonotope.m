function pZ = polyZonotope(SpS)
% polyZonotope - Convert a spectrahedral shadow to a polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    pZ = polyZonotope(SpS);
%
%    figure; hold on;
%    plot(SpS);
%    plot(pZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polyZonotope, taylm/polyZonotope

% Authors:       Adrian Kulmburg
% Written:       06-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pZ = polyZonotope(zonotope(SpS));

% ------------------------------ END OF CODE ------------------------------
