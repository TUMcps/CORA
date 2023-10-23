function example_manual_matZonotope()
% example_manual_matZonotope - example from the manual demonstrating the 
% matZonotope constructor as defined in the manual
%
% Syntax:
%   example_manual_matZonotope()
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

% matrix center
C = [0 0; 0 0];

% matrix generators
G{1} = [1 3; -1 2];
G{2} = [2 0; 1 -1];

% matrix zonotope
mz = matZonotope(C,G);

% ------------------------------ END OF CODE ------------------------------
