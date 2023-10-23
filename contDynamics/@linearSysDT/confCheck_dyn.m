function [res, union_y_a] = confCheck_dyn(sys,params,options)
% confCheck_dyn - specific conformance checker for linear discrete-time
%   systems according to [1] and [2].
%
% Syntax:
%    res = confCheck_dyn(sys,params,options)
%
% Inputs:
%    sys - discrete-time linear system 
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

%% Load variables
nrOfTestCases = length(testSuite); % number of tests
maxNrOfTimeSteps = ceil(params.tFinal/sys.dt); % maximum number of timeSteps

%% Unify all differences to the nominal solution; see y_a(k) in eq. 14 of [1]
[union_y_a,union_x_a] = conform_unifyTestCases(sys,params);

%% Corollary 1 in [2]: change time horizon to ensure conformance of an 
% infinite time horizon; requires that the entire state can be measured
if isinf(params.tFinal)
    % halfspace representation of W and V
    W_halfspace = halfspace(params.W);
    V_halfspace = halfspace(params.V);
    % check enclosure of states after one time step in W
    N = W_halfspace.halfspace.H; % normal vectors
    d = W_halfspace.halfspace.K; % offset of halfspaces
    res_W = all(all(N*union_x_a{2}'<=d)); % k=2 corresponds to t_1
    % check enclosure of deviations of measurements of initial time in V
    N = V_halfspace.halfspace.H; % normal vectors
    d = V_halfspace.halfspace.K; % offset of halfspaces
    % Loop over all test cases; collect differences according to (21) of 
    % [2]; the signs in (21) were incorrect in the IEEE early access version
    union_diff = zeros(nrOfTestCases,size(sys.C,1));
    for iTest=1:nrOfTestCases 
        union_diff(iTest,:) = - sys.C*testSuite{iTest}.x(1,:)' - sys.D*testSuite{iTest}.u(1,:)' + testSuite{iTest}.y(1,:)';
    end
    res_V = all(all(N*union_diff'<=d)); % check enclosure
    % check overall conformance
    res = (res_W && res_V);
else
    % Compute reachable halfspaces, see (eq. (15), (17), Thm. 1 in [1])
    [N, d] = conform_reachableHalfspaces(sys,params,options);

    %% initialize failed test cases
    ind = 1:nrOfTestCases;
    failedCases = cell(1,maxNrOfTimeSteps);
    numFailed = 0;

    %% Test conformance for each time step (Theorem 1 in [1])
    for k = 1:maxNrOfTimeSteps+1
        % get indices of failed test cases through checking halfspace inclusion
        failed = ind(any(~(N{k}*union_y_a{k}'<=d{k}))); % compute indices of failed cases
        numFailed = numFailed + length(failed); % update number of failed cases
        failedCases{k} = failed;
    end
    % conformance is ensured when no measurement is outside the output
    % reachable set
    res = ~numFailed;
end

% ------------------------------ END OF CODE ------------------------------
