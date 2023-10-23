function example_manual_matrixSet_vertices()
% example_manual_matrixSet_vertices - example from the manual demonstrating 
% the vertices operation of a matrix set as defined in the manual
%
% Syntax:
%   example_manual_matrixSet_vertices()
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

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% matrix set
C = [0 1;3 2];
G{1} = [1 2;0 1];
A = matZonotope(C,G);

% compute vertices
res = vertices(A)

% ------------------------------ END OF CODE ------------------------------
