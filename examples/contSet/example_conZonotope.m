function completed = example_conZonotope()
% example_conZonotope - example instantiation of a conZonotope object
%
% Syntax:
%    completed = example_conZonotope()
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

% Authors:       Niklas Kochdumper
% Written:       15-July-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Z = [0 1 0 1; 0 1 2 -1]; % zonotope (center + generators)
A = [-2 1 -1]; % constraints (matrix A)
b = 2; % constraints (vector b)

cZ = conZonotope(Z,A,b) % construct conZonotope object

figure; hold on; box on;
plot(cZ); % visualize conZonotope object
plot(zonotope(cZ),[1,2],'Color',colorblind('r')); % visualize unconstrained zonotope object

completed = true;

% ------------------------------ END OF CODE ------------------------------
