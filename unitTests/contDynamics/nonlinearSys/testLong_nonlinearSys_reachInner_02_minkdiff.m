function res = testLong_nonlinearSys_reachInner_02_minkdiff
% testLong_nonlinearSys_reachInner_02_minkdiff - example for the
%    computation of an inner approximation of the reachable set for
%    nonlinear dynamics, using the Minkowski difference approach from [1]
%
% Syntax:
%    res = testLong_nonlinearSys_reachInner_02_minkdiff
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Wetzlinger, A. Kulmburg, and M. Althoff. "Inner approximations
%        of reachable sets for nonlinear systems using the Minkowski
%        difference.", IEEE Control Systems Letters, 2024.

% Authors:       Mark Wetzlinger
% Written:       17-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% System Dynamics ---------------------------------------------------------

sys = nonlinearSys('jetEngine',@jetEngine); 

% Parameters --------------------------------------------------------------

params.tFinal = 4;
R0 = interval([0.9;0.9],[1.1;1.1]);

% Reachability Settings ---------------------------------------------------

options_outer.alg = 'lin-adaptive';
options_inner.algInner = 'minkdiff';
options_inner.timeStep = 0.01;
options_inner.tensorOrder = 2;
options_inner.compOutputSet = false;

% Reachability Analysis ---------------------------------------------------

params.R0 = zonotope(R0);
Rout = reach(sys,params,options_outer);

params.R0 = polytope(R0);
Rin = reachInner(sys,params,options_inner);


% Validation --------------------------------------------------------------

% check if last inner approximation is contained in outer approximation
Rout_final = query(Rout,'finalSet');
Rin_final = query(Rin,'finalSet');

assert(contains(Rout_final, Rin_final));

% ------------------------------ END OF CODE ------------------------------
