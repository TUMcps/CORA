function example_manual_taylm_composition()
% example_manual_taylm_composition - example from the manual demontrating 
% the composition of Taylor models as defined in the manual
%
% Syntax:
%   example_manual_taylm_composition()
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

a1 = interval(-1, 2); % generate a scalar interval [-1,2]
b1 = taylm(a1, 6); % generate a scalar Taylor model of order 6
a2 = interval(2, 3); % generate a scalar interval [2,3]
b2 = taylm(a2, 6); % generate a scalar Taylor model of order 6
c = [b1; b2] % generate a row of Taylor models

% ------------------------------ END OF CODE ------------------------------
