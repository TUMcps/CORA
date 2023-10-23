function res = testLong_stl_modelCheckReachSet_oscillator()
% testLong_stl_modelCheckReachSet_oscillator - unit test function of
% modelCheckReachSet with an oscillating system
%
% Syntax:
%    res = testLong_stl_modelCheckReachSet_oscillator
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
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

alg = {'signals'}; % 'rtl' and 'sampledTime' are too slow

A = [-0.05 -7; 7 -0.05];
B = 1;
sys = linearSys(A,B);

params.tFinal = 10;
params.R0 = zonotope([[10; 5],0.5*eye(2)]);
params.U = zonotope([ones(2,1),0.25*eye(2)]);

options.linAlg = 'adaptive';
options.error = 0.1;

R = reach(sys,params,options);

x = stl('x',2);

res = true;

% some formulas that must hold

pos{1} = globally(implies(x(1) < -5, finally(x(1) > 5, ...
    interval(0,1))), interval(0,7));

pos{2} = finally(globally(x(1) < 10, interval(0,0.5)),interval(0,9));

for i = 1:length(pos)
    for j = 1:length(alg)
        res = res && modelChecking(R,pos{i},alg{j});
    end
end

% some formulas that must not hold

neg{1} = globally(x(1) > x(2), interval(0,10));
neg{2} = finally(x(1) > 10, interval(1,10));
neg{3} = finally(globally(x(1) > 0, interval(0,1)), interval(0,9));

for i = 1:length(neg)
    for j = 1:length(alg)
        res = res && ~modelChecking(R,neg{i},alg{j});
    end
end

end

% ------------------------------ END OF CODE ------------------------------
