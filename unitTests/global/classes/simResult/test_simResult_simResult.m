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

% init result
res = true;

% empty simResult
try
    simRes = simResult();
catch
    res = false; return
end

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
try
    simRes = simResult(x,t);
    simRes = simResult(x,t,loc);
    simRes = simResult(x,t,{},y);
    simRes = simResult(x,t,{},y,a);
catch
    res = false; return
end

% check wrong instantiations
if CHECKS_ENABLED

% empty time
try
    simRes = simResult(x,{});
    res = false; return
end

% non-matching number of steps
try
    simRes = simResult(x1,t);
    res = false; return
end
try
    simRes = simResult(x,t,{},y1);
    res = false; return
end
try
    simRes = simResult(x,t,{},y,a1);
    res = false; return
end

% non-matching number of runs
try
    simRes = simResult(x2,t,{},y);
    res = false; return
end
try
    simRes = simResult(x,t,{},y2);
    res = false; return
end
try
    simRes = simResult(x,t,{},y,a2);
    res = false; return
end

% different number of columns in between runs
try
    simRes = simResult(x3,t);
    res = false; return
end
try
    simRes = simResult(x,t,{},y3);
    res = false; return
end
try
    simRes = simResult(x,t,{},y,a3);
    res = false; return
end

% too little/many input arguments
try
    simRes = simResult(x);
    res = false; return
end
try
    simRes = simResult(x,t,{},y,a,x);
    res = false; return
end

end

% ------------------------------ END OF CODE ------------------------------
