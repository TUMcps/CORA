function res = test_testCase_sequentialTestCases()
% test_testCase_sequentialTestCases() - unit test function for creating
% sequential test cases of a given length, which start 
% at every instant of the original test case (see Def. 7 in [1])
%
% Syntax:
%    test_testCase_sequentialTestCases()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022

% Authors:       Matthias Althoff
% Written:       15-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'autonomousDriving'];

% load test suite
load([path filesep 'ACC2012Test'],'ACC2012Test');

% create sequential test cases; only for first test case
seqTestCases = sequentialTestCases(ACC2012Test{1}, 20);

% init
resPartial = [];

% check correctness of each test case
for iCase = 1:(length(seqTestCases)-1)
    %% check if values starting from time step 2 equal first values of next test case
    % y
    resPartial(end+1) = all(all(seqTestCases{iCase}.y(2:end,:) == seqTestCases{iCase+1}.y(1:end-1,:)));
    % u
    resPartial(end+1) = all(all(seqTestCases{iCase}.u(2:end,:) == seqTestCases{iCase+1}.u(1:end-1,:)));
    % x
    resPartial(end+1) = all(all(seqTestCases{iCase}.x(2:end,:) == seqTestCases{iCase+1}.x(1:end-1,:)));
    % initial state
    resPartial(end+1) = all(seqTestCases{iCase}.x(2,:) == seqTestCases{iCase+1}.initialState');
end

% final result
res = all(resPartial);
        
% ------------------------------ END OF CODE ------------------------------
