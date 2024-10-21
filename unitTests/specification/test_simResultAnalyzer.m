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

res = true;

% atomic propositions
aps = containers.Map;

aps('p1') = atomicProposition(polytope([0 1], 7));
aps('p2') = atomicProposition(polytope([1 1], 9));

% simulation result
x = [1.1 1; 2.2 2; 3.3 3; 4.4 4; 5.5 5; 6.6 6; 7.7 7; 8.8 8; 9.9 9];
t = [0; 1; 2; 3; 4; 5; 6; 7; 8];
simRes = simResult({x}, {t});

% analyze
a1 = simResultAnalyzer(aps, simRes);

s1 = a1.analyze(1,1);

% test
assert(isequal([7 8], s1('p1').time));
assert(isequal([true false], s1('p1').value));

assert(isequal([4 8], s1('p2').time));
assert(isequal([true false], s1('p2').value));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
