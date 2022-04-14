function completed = example_conZonotope()
% example_conZonotope - example instantiation of a conZonotope objects
%
% Syntax:  
%    completed = example_conZonotope()
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author:        Niklas Kochdumper
% Written:       15-July-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

Z = [0 1 0 1; 0 1 2 -1]; % zonotope (center + generators)
A = [-2 1 -1]; % constraints (matrix A)
b = 2; % constraints (vector b)

cZ = conZonotope(Z,A,b) % construct conZonotope object

plotZono(cZ,[1,2]); % visualize conZonotope object + linear zonotope

completed = 1;

%------------- END OF CODE --------------