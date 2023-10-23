function res = test_reachSetAnalyzer
% test_reachSetAnalyzer - unit test of reachSetAnalyzer class
%
% Syntax:
%    res = test_reachSetAnalyzer
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

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% atomic propositions
aps = containers.Map;

aps('p1') = atomicProposition(halfspace([0 1], 3));
aps('p2') = atomicProposition(halfspace([1 1], 5));


% analyzer
a1 = reachSetAnalyzer(aps, 10);

% reach sets
t1{1} = interval([0.0; 0.0], [1.0; 1.0]);
t1{2} = interval([0.5; 0.5], [1.5; 1.5]);
t1{3} = interval([1.0; 1.0], [2.0; 2.0]);
t1{4} = interval([1.5; 1.5], [2.5; 2.5]);
t1{5} = interval([2.0; 2.0], [3.0; 3.0]);
t1{6} = interval([2.5; 2.5], [3.5; 3.5]);
t1{7} = interval([3.0; 3.0], [4.0; 4.0]);
t1{8} = interval([3.5; 3.5], [4.5; 4.5]);
t1{9} = interval([4.0; 4.0], [5.0; 5.0]);

% analyze
a1.analyzeTimeInterval(t1{1}, 1, interval(0, 2));
a1.analyzeTimeInterval(t1{2}, 1, interval(1, 3));
a1.analyzeTimeInterval(t1{3}, 1, interval(2, 4));
a1.analyzeTimeInterval(t1{4}, 1, interval(3, 5));
a1.analyzeTimeInterval(t1{5}, 1, interval(4, 6));
a1.analyzeTimeInterval(t1{6}, 1, interval(5, 7));
a1.analyzeTimeInterval(t1{7}, 1, interval(6, 8));
a1.analyzeTimeInterval(t1{8}, 1, interval(7, 9));
a1.analyzeTimeInterval(t1{9}, 1, interval(8, 10));

s1 = a1.signals();

% test
res(end+1,1) = isequal([5 8 10], s1('p1').time);
res(end+1,1) = isequal([tt uu ff], s1('p1').value);

res(end+1,1) = isequal([4 7 10], s1('p2').time);
res(end+1,1) = isequal([tt uu ff], s1('p2').value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
