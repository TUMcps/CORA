function example_manual_stl()
% example_manual_stl - example from the manual demontrating the
% stl class as defined in the manual
%
% Syntax:
%   example_manual_stl()
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

% create variable
x = stl('x',2);

% signal temporal logic formula
eq = until(x(1) < 3,x(2) > 5,interval(1,3))

% ------------------------------ END OF CODE ------------------------------
