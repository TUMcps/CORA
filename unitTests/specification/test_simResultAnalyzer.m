function res = test_simResultAnalyzer
% test_simResultAnalyzer - unit test of simResultAnalyzer class
%
% Syntax:
%    res = test_simResultAnalyzer
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Benedikt Seidl
% Written:       24-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% atomic propositions
aps = containers.Map;

aps('p1') = atomicProposition(halfspace([0 1], 7));
aps('p2') = atomicProposition(halfspace([1 1], 9));

% simulation result
x = [1.1 1; 2.2 2; 3.3 3; 4.4 4; 5.5 5; 6.6 6; 7.7 7; 8.8 8; 9.9 9];
t = [0; 1; 2; 3; 4; 5; 6; 7; 8];
simRes = simResult({x}, {t});

% analyze
a1 = simResultAnalyzer(aps, simRes);

s1 = a1.analyze(1);

% test
res(end+1,1) = isequal([7 8], s1('p1').time);
res(end+1,1) = isequal([true false], s1('p1').value);

res(end+1,1) = isequal([4 8], s1('p2').time);
res(end+1,1) = isequal([true false], s1('p2').value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
