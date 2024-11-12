function res = test_linErrorBound_print
% test_linErrorBound_print - unit test for helper class linErrorBound
%
% Syntax:
%    res = test_linErrorBound_print
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

% set partial errors and bounds
errs.timeSteps = {0.1; 0.2; 0.3; 0.3; 0.1};
errs.seq_nonacc = {0.1; 0.15; 0.15; 0.12; 0.05};
errs.bound_rem = {0.5; 0.4; 0.4; 0.2; 0.1};
errs.step_acc = {0.01; 0.01; 0.03; 0.02; 0.02};
errs.bound_acc = {0.1; 0.2; 0.3; 0.3; 0.1};
errs.step_red = [0; 0.04; 0; 0; 0.05];
errs.bound_red = {0.1; 0.2; 0.2; 0.4; 0.5};

% print table
print(errs);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
