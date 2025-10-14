function [res, union_y_a] = isconform_dyn(linsysDT,params,options)
% isconform_dyn - specific conformance checker for linear discrete-time
%   systems according to [1] and [2].
%
% Syntax:
%    res = isconform_dyn(linsysDT,params,options)
%
% Inputs:
%    linsysDT - linearSysDT object
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    res  - true if conformance was achieved, otherwise false
%    union_y_a - unified test cases (only deviation to nominal solution)
%
% References:
%    [1] M. Althoff. Checking and Establishing Reachset Conformance in
%        CORA 2023. In Proc. of the 10th International Workshop on 
%        Applied Verification of Continuous and Hybrid Systems, 2023.
%    [2] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Stefan Liu
% Written:       29-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test suite
testSuite = params.testSuite;

% pre-compute mapping by E and F
params.W = linsysDT.E * params.W;
params.V = linsysDT.F * params.V;

%% Load variables
n_m = length(testSuite); % number of tests
n_k = ceil(params.tFinal/linsysDT.dt); % maximum number of timeSteps

%% Unify all differences to the nominal solution; see y_a(k) in eq. 14 of [1]
[union_y_a,union_x_a] = priv_conform_unifyTrajectories(linsysDT,params);

%% Corollary 1 in [2]: change time horizon to ensure conformance of an 
% infinite time horizon; requires that the entire state can be measured
if isinf(params.tFinal)
    % halfspace representation of W and V
    W_polytope = polytope(params.W);
    V_polytope = polytope(params.V);
    % check enclosure of states after one time step in W
    N = W_polytope.A; % normal vectors
    d = W_polytope.b; % offset of halfspaces
    res_W = all(all(N*union_x_a{2}<=d)); % k=2 corresponds to t_1

    % check enclosure of deviations of measurements of initial time in V
    N = V_polytope.A; % normal vectors
    d = V_polytope.b; % offset of halfspaces
    % Loop over all test cases; collect differences according to (21) of 
    % [2]; the signs in (21) were incorrect in the IEEE early access version
    union_diff = zeros(n_m,size(linsysDT.C,1));
    for m=1:n_m 
        union_diff(m,:) = - linsysDT.C*testSuite(m).x(:,1) - linsysDT.D*testSuite(m).u(:,1) + testSuite(m).y(:,1);
    end
    res_V = all(all(N*union_diff'<=d)); % check enclosure
    
    % check overall conformance
    res = (res_W && res_V);
else
    % Compute reachable halfspaces, see (eq. (15), (17), Thm. 1 in [1])
    [N, d] = priv_conform_reachableHalfspaces(linsysDT,params,options);

    %% initialize failed test cases
    ind = 1:n_m*testSuite(1).n_s;
    failedCases = cell(1,n_k);
    numFailed = 0;

    %% Test conformance for each time step (Theorem 1 in [1])
    for k = 1:n_k+1
        % get indices of failed test cases through checking halfspace inclusion
        failed = ind(any(~(N{k}*union_y_a{k}<=d{k}))); % compute indices of failed cases
        numFailed = numFailed + length(failed); % update number of failed cases
        failedCases{k} = failed;
    end
    % conformance is ensured when no measurement is outside the output
    % reachable set
    res = ~numFailed;
end

% ------------------------------ END OF CODE ------------------------------
