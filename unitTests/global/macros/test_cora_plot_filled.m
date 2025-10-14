function res = test_cora_plot_filled
% test_cora_plot_filled - tests the macro CORA_PLOT_FILLED
%
% Syntax:
%    res = test_cora_plot_filled()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CHECKS_ENABLED

% Authors:       Tobias Ladner
% Written:       11-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

assert(CORA_PLOT_FILLED);

% ------------------------------ END OF CODE ------------------------------
