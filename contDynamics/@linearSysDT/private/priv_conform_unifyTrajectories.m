function [union_y_a, union_x_a] = priv_conform_unifyTrajectories(linsysDT,params)
% priv_conform_unifyTrajectories - unifies trajectories according to 
%    Sec. III.B in [1] and Sec. IV.C in [2].
%
% Syntax:
%    res = priv_conform_unifyTrajectories(linsysDT,params,options)
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
%                08-September-2025 (LL, adapt to trajectories)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test suite
testSuite = params.testSuite;

%% Load variables
n_m = length(testSuite); % number of tests
n_s = size(testSuite(1).y,3);
dim_y = size(linsysDT.C,1); % output dimension
dim_u = size(linsysDT.B,2); % input dimension
dim_x = size(linsysDT.A,1); % state dimension

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
[union_y_a{:}] = deal(zeros(dim_y,n_m*n_s));
[union_x_a{:}] = deal(zeros(dim_x,n_m*n_s));

% Loop over all test cases
for m=1:n_m
    % check if time step size is correct
    assert(abs(linsysDT.dt-testSuite(m).dt)<1e-6,"Test case number "+m+" has an incorrect sampling time ("+testSuite(m).dt+"). Correct sampling rate should be "+linsysDT.dt);

    for s=1:n_s

        % set initial state and final time of current test case
        n_k = testSuite(m).n_k; % timeSteps
        params.tFinal = (n_k-1)*linsysDT.dt;
        params.x0 = testSuite(m).x(:,1,1);

        %% set input trajectory
        % no input
        if isempty(testSuite(m).u)
            params.u = zeros(dim_u);
            % input exists
        else
            params.u = testSuite(m).u(:,1:n_k);
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
        y_a = testSuite(m).y(:,1:n_k,s) - y_star;
        % Compute difference to nominal state solution (if it is provided)
        try
            x_a = testSuite(m).x(:,1:n_k,s) - x_star;
        catch
            x_a = [];
        end
        % restructure output differences so that differences at the same time
        % are collected
        for k = 1:n_k
            union_y_a{k}(:,(m-1)*n_s+s) = y_a(:,k);
        end
        % restructure state differences so that differences at the same time
        % are collected
        if ~isempty(x_a)
            for k = 1:n_k
                union_x_a{k}(:,(m-1)*n_s+s) = x_a(:,k);
            end
        else
            union_x_a(:) = [];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
