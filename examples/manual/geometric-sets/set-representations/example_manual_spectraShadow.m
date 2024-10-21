function example_manual_spectraShadow()
% example_manual_spectraShadow - example from the manual demonstrating the 
% spectraShadow constructor as defined in the manual
%
% Syntax:
%   example_manual_spectraShadow()
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

enlargeAxis(1.1)
title('$\mathcal{SPS}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
