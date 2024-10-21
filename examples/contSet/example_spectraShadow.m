function completed = example_spectraShadow()
% example_spectraShadow - example instantiation of spectraShadow objects
%
% Syntax:
%    completed = example_spectraShadow()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       08-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create two simple spectrahedral shadows
% an ellipsoid, using the first type of instantiation
A0 = eye(3);
A1 = [0 1 0; 1 0 0; 0 0 0];
A2 = [0 0 1; 0 0 0; 1 0 0];
A = [A0 A1 A2];
c = [-1.5;0];
G = [1 0;0 1.5];

SpS_ellipsoid = spectraShadow(A,c,G);

% a small zonotope, using the second type of instantiation
A0 = diag([-19 1 11 21 1 -9]);
A1 = diag([10 10 0 -10 -10 0]);
A2 = (10/3)*diag([1 -1 2 -1 1 -2]);
A = [A0 A1 A2];

B = 0.5774*diag([-1 1 1 1 -1 -1]);

ESumRep = {A, B};

SpS_zonotope = spectraShadow(ESumRep);

% we can construct more complicated sets
% based on the two above, such as the
% convex hull of the two previous spectrahedral shadows
SpS = convHull(SpS_ellipsoid, SpS_zonotope);

figure;
plot(SpS)

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
SpS_addition = SpS_polytope + SpS_ellipsoid 
% compute convex hull of C and Z
SpS_convHull = convHull(SpS_capsule, SpS_zonotope); 

figure; hold on;
% One can now plot the convex hull of C and Z
plot(SpS_convHull,[1 2],'FaceColor',colorblind('gray'));
plot(C,[1 2],'Color',colorblind('r')); % plot C
plot(Z,[1 2],'Color',colorblind('b')); % plot Z


%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
