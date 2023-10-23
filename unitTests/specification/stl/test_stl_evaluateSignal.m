function res = test_stl_evaluateSignal
% test_stl_evaluateSignal - unit test of stl evaluateSignal method
%
% Syntax:
%    res = test_stl_evaluateSignal
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

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% atomic propositions
p = stl('P', atomicProposition(interval(0,1)));
q = stl('Q', atomicProposition(interval(0,1)));

sigs = containers.Map;

% signals
t1 = [1.5 2.5 7.0 8.5 10];
v1 = [ff uu tt uu ff];
sigs("P") = signal(t1, v1);

t2 = [2.0 3.5 4.5 5.5 6.5 7.5 9.0 10];
v2 = [ff uu tt uu tt ff uu ff];
sigs("Q") = signal(t2, v2);

% test
s1 = evaluateSignal(until(p, q, interval(0,0.5)), 9, sigs);
res = isequal([1.5 3.0 4.5 5.0 6.5 7.0 8.5 9], s1.time);
res(end+1,1) = isequal([ff uu tt uu tt ff uu ff], s1.value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
