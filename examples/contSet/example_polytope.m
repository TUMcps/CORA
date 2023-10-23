function completed = example_polytope()
% example_polytope - example instantiation of polytope objects
%
% Syntax:
%    completed = example_polytope()
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

% Authors:       Matthias Althoff
% Written:       21-April-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([-1 1 0; 1 0 1]); % create zonotope Z2

P1 = polytope(Z1); % convert zonotope Z1 to halfspace representation
P2 = polytope(Z2); % convert zonotope Z2 to halfspace representation

P3 = P1 + P2 % perform Minkowski addition and display result
P4 = P1 & P2; % compute intersection of P1 and P2

V = vertices(P4) % obtain and display vertices of P4

figure; hold on;
plot(P4,[1 2],'FaceColor',colorblind('gray')); % plot P4 
plot(P3,[1 2],'Color',colorblind('r')); % plot P3
plot(P1); % plot P1
plot(P2); % plot P2


%example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
