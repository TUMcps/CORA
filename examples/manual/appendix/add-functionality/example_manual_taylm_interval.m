function example_manual_taylm_interval()
% example_manual_taylm_interval - example from the manual demontrating 
% the conversion from an interval to a Taylor models as defined in the manual
%
% Syntax:
%   example_manual_taylm_interval()
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

a = interval([-1;2], [2;3]); % generate an interval vector [[-1,2]; [2,3]]
c = taylm(a, 6, {'a1';'a2'}) % generate Taylor model (order 6)

% ------------------------------ END OF CODE ------------------------------
