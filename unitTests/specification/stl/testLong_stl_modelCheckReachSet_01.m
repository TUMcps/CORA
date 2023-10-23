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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
res = true;
alg = {'rtl','sampledTime','signals'};

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
for j = 1:length(alg)
    for i = 1:4
    
        % reachability analysis
        options.timeStep = options.timeStep/2;
    
        R = reach(sys,params,options);
    
        % model checking (should be false for arbitrary time steps)
        if modelChecking(R,eq,alg{j}) ~= false
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% modified STL formula
y = stl('y',2);
point = 0.7;
eq = until(y(2) < -point,y(1) > point,interval(0,2));

% refine time step until equation satisfied
options.timeStep = 1;

for j = 1:length(alg)

    resTmp = false;

    for i = 1:5
    
        % reachability analysis
        options.timeStep = options.timeStep/2;
    
        R = reach(sys,params,options);
    
        % model checking (should become true for small enought time steps)
        if modelChecking(R,eq,alg{j}) == true
            resTmp = true; break;
        end
    end

    if ~resTmp
        throw(CORAerror('CORA:testFailed'));
    end
end

% ------------------------------ END OF CODE ------------------------------
