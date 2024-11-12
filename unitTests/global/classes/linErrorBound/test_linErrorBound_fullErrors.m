function res = test_linErrorBound_fullErrors
% test_linErrorBound_fullErrors - unit test for helper class linErrorBound
%
% Syntax:
%    res = test_linErrorBound_fullErrors
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

% tolerance
tol = 1e-14;

% init
emax = 10; tFinal = 1;
errs = linErrorBound(emax,tFinal);

% set non-accumulating/accumulating/reduction errors
errs.seq_nonacc = {1;2;3;4;5};
errs.cum_acc = 0.1*(1:5)';
errs.cum_red = 0.01*(1:5)';

% compute errors of time-point and time-interval solution
[Rcont_error,Rcont_tp_error] = fullErrors(errs,1);
assert(withinTol(Rcont_error,1.01,tol));
assert(withinTol(Rcont_tp_error,0.11,tol));

[Rcont_error,Rcont_tp_error] = fullErrors(errs,5);
assert(withinTol(Rcont_error,5.45,tol));
assert(withinTol(Rcont_tp_error,0.55,tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
