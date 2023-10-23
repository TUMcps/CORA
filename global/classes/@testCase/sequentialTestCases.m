function res = sequentialTestCases(obj, length)
% sequentialTestCases - Creates sequential test cases of a given length, 
% starting at every instant of the original test case (see Def. 7 in [1])
%
% Syntax:
%    obj = sequentialTestCases(obj,length)
%
% Inputs:
%    obj - testCase object
%    length - number of samples for each test case
%
% Outputs:
%    res - cell array of testCase objects
%
% References:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022
%
% Example:
%    ACC2012Test = load(testCase_ACC2012);
%    res = sequentialTestCases(ACC2012Test{1}, 20)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff, Stefan Liu
% Written:       15-June-2023             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% sanity checks
% check if testcase object is empty
if isempty(obj)
    res{1} = obj;
    return
end
% check if states are available
assert(~isempty(obj.x),"To create sequential test cases, the state must be available at every sampling instant. Please add the states to the test case before running this method.");
% check if length of test case is long enough
nrOfTestCases = size(obj.y,1)-length;
if nrOfTestCases < 1
    res{1} = obj;
    return
end

% create template object for each new test case
dataStruct = obj;

% initialize result
res = cell(nrOfTestCases,1);

% loop for each new test case
for i = 1:nrOfTestCases
    % update measured outputs
    dataStruct.y = obj.y(i:i+length,:);
    % update inputs
    if ~isempty(obj.u)
        dataStruct.u = obj.u(i:i+length,:);
    end
    % update states
    dataStruct.x = obj.x(i:i+length,:);
    % update initial state
    dataStruct.initialState = obj.x(i,:)';
    % save result in cell array
    res{i} = dataStruct;
end
end

% ------------------------------ END OF CODE ------------------------------
