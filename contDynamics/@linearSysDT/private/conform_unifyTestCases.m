function [union_y_a, union_x_a] = conform_unifyTestCases(linsysDT,params)
% conform_unifyTestCases - unifies test cases according to Sec. III.B in 
%    [1] and Sec. IV.C in [2].
%
% Syntax:
%    res = conform_unifyTestCases(linsysDT,params,options)
%
% Inputs:
%    linsysDT - linearSysDT object
%    params - parameter defining the conformance problem
%
% Outputs:
%    union_y_a - union all output differences to the nominal solution; see 
%                y_a(k) in eq. 14 of [1]
%    union_x_a - union all state differences to the nominal solution; does
%                not exist if states are not contained in test cases
%
% Example: 
%    -
%
% References:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%    [2] S. B. Liu and M. Althoff. Reachset conformance of forward dynamic 
%        models for the formal analysis of robots. In Proc. of the IEEE/RSJ 
%        International Conference on Intelligent Robots and Systems, 
%        page 370–376, 2018.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-June-2023
% Last update:   14-May-2024 (LL, adapt to 3D arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test suite
testSuite = params.testSuite;

%% Load variables
nrOfTestCases = length(testSuite); % number of tests
nrOfSamples = size(testSuite{1}.y,3);
q = size(linsysDT.C,1); % output dimension
p = size(linsysDT.B,2); % input dimension
n = size(linsysDT.A,1); % state dimension

%% Unify all differences to the nominal solution; see y_a(k) in eq. 14 of [1]
% initialize union of differences to the nominal solution
% if the time horizon is infinite, we do not fully initialize all variables
if isinf(params.tFinal)
    maxNrOfTimeSteps = 1;
else
    maxNrOfTimeSteps = ceil(params.tFinal/linsysDT.dt); % maximum number of timeSteps
end
union_y_a = cell(maxNrOfTimeSteps,1);
union_x_a = cell(maxNrOfTimeSteps,1);
[union_y_a{:}] = deal(zeros(nrOfTestCases*nrOfSamples,q));
[union_x_a{:}] = deal(zeros(nrOfTestCases*nrOfSamples,n));

% Loop over all test cases
for iTest=1:nrOfTestCases
    % check if time step size is correct
    assert(linsysDT.dt==testSuite{iTest}.sampleTime,"Test case number "+iTest+" has an incorrect sampling time ("+testSuite{iTest}.sampleTime+"). Correct sampling rate should be "+linsysDT.dt);

    for is=1:nrOfSamples

        % set initial state and final time of current test case
        nrOfTimeSteps = size(testSuite{iTest}.y(:,:,is),1); % timeSteps
        params.tFinal = (nrOfTimeSteps-1)*linsysDT.dt;
        params.x0 = testSuite{iTest}.initialState(:,:,is)';

        %% set input trajectory
        % no input
        if isempty(testSuite{iTest}.u)
            params.u = zeros(p);
            % input exists
        else
            params.u = testSuite{iTest}.u(1:nrOfTimeSteps,:,is)';
            % the length of the input trajectory has to be reduced by one if there
            % is no D to properly compute the output, see simulation function of
            % this class
            if ~any(any(linsysDT.D))
                params.u(:,end) = [];
            end
        end

        %% compute difference to nominal solution
        % nominal solution
        [~,x_star,~,y_star] = simulate(linsysDT,params);
        % Compute difference to nominal output solution
        y_a = testSuite{iTest}.y(1:nrOfTimeSteps,:,is) - y_star;
        % Compute difference to nominal state solution (if it is provided)
        try
            x_a = testSuite{iTest}.x(1:nrOfTimeSteps,:,is) - x_star;
        catch
            x_a = [];
        end
        % restructure output differences so that differences at the same time
        % are collected
        for k = 1:nrOfTimeSteps
            union_y_a{k}((iTest-1)*nrOfSamples+is,:) = y_a(k,:);
        end
        % restructure state differences so that differences at the same time
        % are collected
        if ~isempty(x_a)
            for k = 1:nrOfTimeSteps
                union_x_a{k}((iTest-1)*nrOfSamples+is,:) = x_a(k,:);
            end
        else
            union_x_a(:) = [];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
