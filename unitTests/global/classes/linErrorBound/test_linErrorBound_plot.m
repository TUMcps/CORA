function res = test_linErrorBound_plot
% test_linErrorBound_plot - unit test for helper class linErrorBound
%
% Syntax:
%    res = test_linErrorBound_plot
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
emax = 10; tFinal = 10;
errs = linErrorBound(emax,tFinal);

% set partial errors and bounds
errs.timeSteps = {1;2;1;3;2;1};
errs.seq_nonacc = {4;3.5;2;2.75;2;1};
errs.bound_rem = {5;4;3.5;3;2.75;2.5};
errs.step_acc = {0.1;0.3;0.2;0.5;0.3;0.1};
errs.bound_acc = {0.2;0.5;0.3;1;0.5;0.2};
errs.step_red = [0;0;0.2;0.5;0;0.1];
errs.bound_red = {0.5;1.5;2;3;2;1};

% plot errors
plot(errs);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
