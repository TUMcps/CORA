function res = testLong_nonlinearSys_reach_03_vanDerPol
% testLong_nonlinearSys_reach_03_vanDerPol - unit_test_function of nonlinear
%    reachability analysis; Checks the solution of the nonlinearSys class
%    for the van der Pol example; The settings are identical to [1].
%    It is checked whether the reachable set is enclosed in the initial set
%    after a certain amount of time.
%
% Syntax:
%    res = testLong_nonlinearSys_reach_03_vanDerPol
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Althoff, O. Stursberg, M. Buss
%        "Reachability analysis of nonlinear systems with uncertain
%        parameters using conservative linearization", CDC 2008

% Authors:       Matthias Althoff
% Written:       26-June-2009
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 6.74;
Z0{1} = zonotope([1.4; 2.3],[0.05 0; 0 0.05]);
params.R0 = zonoBundle(Z0);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.02;
options.taylorTerms = 4;
options.zonotopeOrder = 10;
options.intermediateOrder = 10;
options.errorOrder = 5;

options.alg = 'lin';
options.tensorOrder = 3;

options.maxError = 0.05*[1; 1];
options.reductionInterval = 100;


% System Dynamics ---------------------------------------------------------

vanderPol=nonlinearSys(@vanderPolEq); %initialize van-der-Pol oscillator


% Reachability Analysis ---------------------------------------------------

R = reach(vanderPol, params, options);


% Verification ------------------------------------------------------------

% obtain array of enclosing polytopes of last reachable set
Rfin = query(R,'finalSet');
Premain = polytope(Rfin{1});
for i = 2:length(Rfin)
    Premain = Premain | polytope(Rfin{i});
end

% check containment in hull of first couple of time-interval solutions
P_hull = polytope(params.R0);
for s=1:6
    P_hull = P_hull | polytope(R(1).timeInterval.set{s});
end
res = contains(P_hull,Premain);


% fix polytope/mldivide for test below...
% % remove previous reachable sets
% iStep = 1; t = 0;
% while ~isempty(Premain) && t < 0.4
%     for iSet = 1:length(R(iStep).timeInterval.set)
%         Preach = polytope(R(iStep).timeInterval.set{iSet});
%         Premain = Premain \ Preach;
%         t = R(iStep).timePoint.time{iSet};
%         if t > 0.4
%            break; 
%         end
%     end
%     iStep = iStep + 1;
% end
% 
% % obtain result
% if representsa(Premain,'emptySet')
%     res = true;
% else
%     res = false;
% end

% ------------------------------ END OF CODE ------------------------------
