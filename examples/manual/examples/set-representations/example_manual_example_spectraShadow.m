function example_manual_example_spectraShadow()
% example_manual_example_spectraShadow - example from the manual 
% demonstrating the spectrahedral shadow example from the manual
%
% Syntax:
%   example_manual_example_spectraShadow()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Adrian Kulmburg
% Written:       10-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Any convex set representation implemented in CORA can be represented as a
% spectrahedral shadow. For instance, ellipsoids, capsules, polytopes, and
% zonotopes can be recast as spectrahedral shadows:
E = ellipsoid([2 1; 1 2], [3;1]);
SpS_ellipsoid = spectraShadow(E);

C = capsule([-1;-3], [1;-2], 1);
SpS_capsule = spectraShadow(C);

P = polytope([1 0 -1 0 1; 0 1 0 -1 1]', [3; 2; 3; 2; 1]);
SpS_polytope = spectraShadow(P);

Z = zonotope([3;2], [0.5 1 0; 0 0.5 2.5]);
SpS_zonotope = spectraShadow(Z);

% perform Minkowski addition and display result
SpS_addition = SpS_polytope + SpS_ellipsoid; 
% compute convex hull of C and Z
SpS_convHull = convHull(SpS_capsule, SpS_zonotope); 

figure; hold on;
% One can now plot the convex hull of C and Z
plot(SpS_convHull,[1 2],'FaceColor',colorblind('gray'));
plot(C,[1 2],'Color',colorblind('r')); % plot C
plot(Z,[1 2],'Color',colorblind('b')); % plot Z

% plot --------------------------------------------------------------------

enlargeAxis(1.2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
