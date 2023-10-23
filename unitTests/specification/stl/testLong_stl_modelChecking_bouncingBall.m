function res = testLong_stl_modelChecking_bouncingBall
% testLong_stl_modelChecking_bouncingBall - unit test function of
%    modelChecking
%
% Syntax:
%    res = testLong_stl_modelChecking_bouncingBall
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

% Authors:       Benedikt Seidl
% Written:       15-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

alg = {'rtl','sampledTime','signals'};

% hybrid system
sys = bouncing_ball(-0.75);

params.R0 = zonotope([1;0],diag([0.05,0.05]));
params.startLoc = 1;
params.tFinal = 2;

% continuous dynamics
options.timeStep = 0.01;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% hybrid dynamics
options.guardIntersect = 'polytope';
options.enclose = {'box'};

R = reach(sys, params, options);

% STL formulas
x = stl('x',2);

% loop over different formulas and algorithms
res = [];


% formulas to verify
eq_pos = {
    finally(x(1) > 0.5, interval(0.5,1))
    finally(x(1) < 0.1, interval(0,1))
    finally(x(1) < 0.5, interval(0,1)) & ...
        globally(x(1) < 1.1, interval(0,1))
    globally(finally(x(2) < 0, interval(0,1)), interval(0,1))
    globally(implies(x(1) < 0.1, finally(x(1) > 0.1, ...
        interval(0,0.2))), interval(0,1))
    globally(x(1) < 1.1, interval(0,2))
    finally(x(1) < 0.5, interval(0,1))
    finally(x(1) < 0.1, interval(0,2))
    globally(implies(x(1) < 0.1, finally(x(1) > 0.1, ...
        interval(0,0.5))), interval(0,1.5))
};

for j = 3:length(alg)
    for i = 1:length(eq_pos)
        res(end+1,1) = modelChecking(R,eq_pos{i},alg{j});
    end
end


% formulas to falsify
eq_neg = {
    globally(x(1) < 0.1, interval(0,1))
    finally(x(1) < -0.1, interval(0,1))
    globally(finally(x(1) > 0.4, interval(0,1)), interval(0,1))
    globally(x(1) < 0, interval(0,1))
    globally(x(1) < 0.1, interval(0,1))
    globally(x(1) > 0.1, interval(0,1))
    finally(x(1) > 1.1, interval(0,2))
    finally(x(1) < -0.1, interval(0,1))
    finally(x(1) < 0.5, interval(0,0.2))
    finally(x(1) > 2, interval(0,1))
    globally(finally(x(1) > 0.4, interval(0,1)), interval(0,1))
};

for j = 3:length(alg)
    for i = 1:length(eq_neg)
        res(end+1,1) = ~modelChecking(R,eq_neg{i},alg{j});
    end
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
