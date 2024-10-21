function res = test_simResult_simResult
% test_simResult_simResult - unit test function for constructor
%
% Syntax:
%    res = test_simResult_simResult()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       11-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty simResult
simRes = simResult();

% initialize some sets for instantiations
steps = 10;
runs = 5;
t = {cumsum([0;rand(steps,1)])};

nrStates = 4; nrOutputs = 2; nrAlgebraic = 3;
x = cell(runs,1); x1 = cell(runs,1); x2 = cell(runs-1,1); x3 = cell(runs,1);
y = cell(runs,1); y1 = cell(runs,1); y2 = cell(runs-1,1); y3 = cell(runs,1);
a = cell(runs,1); a1 = cell(runs,1); a2 = cell(runs-1,1); a3 = cell(runs,1);
loc = {1};

% instantiate trajectories
for i=1:runs
    x{i,1} = randn(steps+1,nrStates);
    y{i,1} = randn(steps+1,nrOutputs);
    a{i,1} = randn(steps+1,nrAlgebraic);
    % wrong number of steps
    x1{i,1} = randn(steps,nrStates);
    y1{i,1} = randn(steps,nrOutputs);
    a1{i,1} = randn(steps,nrAlgebraic);
    % wrong number of runs
    if i < runs
        x2{i,1} = randn(steps+1,nrStates);
        y2{i,1} = randn(steps+1,nrOutputs);
        a2{i,1} = randn(steps+1,nrAlgebraic);
    end
    % different number of states/outputs/algebraics in between runs
    if i == 2
        x3{i,1} = randn(steps,nrStates+1);
        y3{i,1} = randn(steps,nrOutputs+1);
        a3{i,1} = randn(steps,nrAlgebraic+1);
    else
        x3{i,1} = randn(steps,nrStates);
        y3{i,1} = randn(steps,nrOutputs);
        a3{i,1} = randn(steps,nrAlgebraic);
    end
end

% correct instantiations according to constructor (all with correct length)
simRes = simResult(x,t);
simRes = simResult(x,t,loc);
simRes = simResult(x,t,{},y);
simRes = simResult(x,t,{},y,a);

% check wrong instantiations

% empty time
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,{});

% non-matching number of steps
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x1,t);
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,t,{},y1);
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,t,{},y,a1);

% non-matching number of runs
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x2,t,{},y);
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,t,{},y2);
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,t,{},y,a2);

% different number of columns in between runs
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x3,t);
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,t,{},y3);
assertThrowsAs(@simResult,'CORA:wrongInputInConstructor',x,t,{},y,a3);

% too many input arguments
assertThrowsAs(@simResult,'CORA:numInputArgsConstructor',x,t,{},y,a,x);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
