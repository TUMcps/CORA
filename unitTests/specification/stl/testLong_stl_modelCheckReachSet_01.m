function res = testLong_stl_modelCheckReachSet_01()
% testLong_stl_modelCheckReachSet_01 - unit test function of
%    modelCheckReachSet
%
% Syntax:
%    res = testLong_stl_modelCheckReachSet_01
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
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   21-February-2024 (FL, add incremental algorithm and speed up by reusing reach set)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
res = true;
alg = {'rtl','sampledTime','signals','incremental'};

% analytical test for the until-operator

% dynamic system, parameters, and options
sys = linearSys([0 -1; 1 0],[0;0]);

params.R0 = zonotope([0;-1]);
params.tFinal = 2;

options.timeStep = 1;
options.zonotopeOrder = 10;
options.taylorTerms = 10;

% STL formula
y = stl('y',2);
point = sqrt(2)/2;
eq = until(y(2) < -point,y(1) > point,interval(0,2));

% loop over different time step sizes and algorithms
for i = 1:4
    % reachability analysis
    options.timeStep = options.timeStep/2;
    R = reach(sys,params,options);

    for j = 1:length(alg)
        % model checking (should be false for arbitrary time steps)
        assertLoop(~modelChecking(R,eq,alg{j}),i,j)
    end
end

% modified STL formula
y = stl('y',2);
point = 0.7;
eq = until(y(2) < -point,y(1) > point,interval(0,2));

% refine time step until equation satisfied
options.timeStep = 1;


for i = 1:5
    resTmp = false(length(alg),1);
    % reachability analysis
    options.timeStep = options.timeStep/2;
    R = reach(sys,params,options);

    for j = 1:length(alg)
        % model checking (should become true for small enough time steps)
        resTmp(j) = modelChecking(R,eq,alg{j});
    end

    if all(resTmp)
        break;
    end
end

% if still not true...
assert(all(resTmp))

% ------------------------------ END OF CODE ------------------------------
