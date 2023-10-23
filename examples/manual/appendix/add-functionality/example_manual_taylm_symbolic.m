function example_manual_taylm_symbolic()
% example_manual_taylm_symbolic - example from the manual demontrating 
% the instantiate Taylor models using symbolic expression as defined in the manual
%
% Syntax:
%   example_manual_taylm_symbolic()
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

syms a1 a2; % instantiate symbolic variables
s = [2 + 1.5*a1; 2.75 + 0.25*a2]; % create symbolic function
c = taylm(s, interval([-2;-3],[0;1]), 6) % generate Taylor model

% ------------------------------ END OF CODE ------------------------------
