function zB = zonoBundle(SpS)
% zonoBundle - overapproximates a spectrahedral shadow by a zonotope bundle
%
% Syntax:
%    zB = zonoBundle(SpS)
%
% Inputs:
%    SpS - spectrahedron object
%
% Outputs:
%    zB - zonoBundle object
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    zB = zonoBundle(SpS);
%
%    figure; hold on;
%    plot(SpS,[1,2],'b');
%    plot(zB,[1,2],'r');
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

zB = zonoBundle(zonotope(SpS));

% ------------------------------ END OF CODE ------------------------------
