function example_manual_intervalMatrix()
% example_manual_intervalMatrix - example from the manual demonstrating the 
% intervalMatrix constructor as defined in the manual
%
% Syntax:
%   example_manual_intervalMatrix()
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

% center matrix
C = [0 2; 3 1];

% width matrix
D = [1 2; 1 1];

% interval matrix
mi = intervalMatrix(C,D)

% ------------------------------ END OF CODE ------------------------------
