function Z = zonotope(SpS)
% zonotope - overapproximates a spectrahedral shadow by a zonotope
%
% Syntax:
%    Z = zonotope(SpS)
%
% Inputs:
%    S - spectraShadow object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    Z = zonotope(SpS);
%
%    figure; hold on;
%    plot(SpS,[1,2],'b');
%    plot(Z,[1,2],'r');
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

% cover the empty case
if representsa_(SpS,'emptySet',1e-5)
    Z = zonotope.empty(dim(SpS)); 
    return
end


% extract the 'base' spectrahedron that is not projected
SpS_base = spectraShadow(SpS.A);

% transform the base spectrahedron to a box
I = interval(SpS_base);

% transform that box to a zonotope
Z = zonotope(I);

% and finally apply the affine transformation of S
Z = full(SpS.G) * Z + full(SpS.c);

% ------------------------------ END OF CODE ------------------------------
