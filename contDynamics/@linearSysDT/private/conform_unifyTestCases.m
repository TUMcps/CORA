function [union_y_a, union_x_a] = conform_unifyTestCases(obj,params)
% conform_unifyTestCases - unifies test cases according to Sec. III.B in 
% [1] and Sec. IV.C in [2].
%
% Syntax:
%    res = conform_unifyTestCases(obj,params,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    params - parameter defining the conformance problem
%
% Outputs:
%    union_y_a  - union all output differences to the nominal solution; see 
%                 y_a(k) in eq. 14 of [1]
%    union_x_a  - union all state differences to the nominal solution; does
%                 not exist if states are not contained in test cases
%
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%    [2] S. B. Liu and M. Althoff. Reachset conformance of forward dynamic 
%        models for the formal analysis of robots. In Proc. of the IEEE/RSJ 
%        International Conference on Intelligent Robots and Systems, 
%        page 370–376, 2018.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test suite

testSuite = params.testSuite;

%% Load variables
nrOfTestCases = length(testSuite); % number of tests
q = size(obj.C,1); % output dimension
p = size(obj.B,2); % input dimension
n = size(obj.A,1); % state dimension

%% Unify all differences to the nominal solution; see y_a(k) in eq. 14 of [1]
% initialize union of differences to the nominal solution
% if the time horizon is infinite, we do not fully initialize all variables
if isinf(params.tFinal)
    maxNrOfTimeSteps = 1;
else
    maxNrOfTimeSteps = ceil(params.tFinal/obj.dt); % maximum number of timeSteps
end
union_y_a = cell(maxNrOfTimeSteps,1);
union_x_a = cell(maxNrOfTimeSteps,1);
[union_y_a{:}] = deal(zeros(nrOfTestCases,q));
[union_x_a{:}] = deal(zeros(nrOfTestCases,n));

% Loop over all test cases
for iTest=1:nrOfTestCases 
    
    % check if time step size is correct
    assert(obj.dt==testSuite{iTest}.sampleTime,"Test case number "+iTest+" has an incorrect sampling time ("+testSuite{iTest}.sampleTime+"). Correct sampling rate should be "+obj.dt);
    
    % set initial state and final time of current test case
    nrOfTimeSteps = size(testSuite{iTest}.y,1); % timeSteps
    params.tFinal = (nrOfTimeSteps-1)*obj.dt;
    params.x0 = testSuite{iTest}.initialState';
    
    %% set input trajectory
    % no input
    if isempty(testSuite{iTest}.u)
        params.u = zeros(p);
    % input exists
    else
        params.u = testSuite{iTest}.u(1:nrOfTimeSteps,:)'; 
        % the length of the input trajectory has to be reduced by one if there 
        % is no D to properly compute the output, see simulation function of
        % this class
        if ~any(any(obj.D))
            params.u(:,end) = [];
        end
    end     
    
    %% compute difference to nominal solution
    % nominal solution
    [~,x_star,~,y_star] = simulate(obj,params); 
    % Compute difference to nominal output solution
    y_a = testSuite{iTest}.y(1:nrOfTimeSteps,:) - y_star; 
    % Compute difference to nominal state solution (if it is provided)
    try
        x_a = testSuite{iTest}.x(1:nrOfTimeSteps,:) - x_star; 
    catch
        x_a = [];
    end
    % restructure output differences so that differences at the same time 
    % are collected
    for k = 1:nrOfTimeSteps
        union_y_a{k}(iTest,:) = y_a(k,:); 
    end
    % restructure state differences so that differences at the same time 
    % are collected
    if ~isempty(x_a)
        for k = 1:nrOfTimeSteps
            union_x_a{k}(iTest,:) = x_a(k,:); 
        end
    else
        union_x_a(:) = [];
    end
end


% ------------------------------ END OF CODE ------------------------------
