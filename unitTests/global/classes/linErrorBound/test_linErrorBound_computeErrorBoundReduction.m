function res = test_linErrorBound_computeErrorBoundReduction
% test_linErrorBound_computeErrorBoundReduction - unit test for helper
%    class linErrorBound
%
% Syntax:
%    res = test_linErrorBound_computeErrorBoundReduction
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       10-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
emax = 10; tFinal = 1;
errs = linErrorBound(emax,tFinal);
A = [-1 4; -4 -1];

% compute allocation for reduction error
% no inputs: no reduction error required in this case
G_U = [0; 0];
computeErrorBoundReduction(errs,A,G_U);
assert(errs.bound_red_max == 0);

% ...with inputs: must be larger than 0 and smaller than maximum error
G_U = [1 0; -0.5 0];
computeErrorBoundReduction(errs,A,G_U);
assert(errs.bound_red_max > 0 && errs.bound_red_max < errs.emax);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
