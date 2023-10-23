function res = test_stl_desugar
% test_stl_desugar - unit test of stl desugar method
%
% Syntax:
%    res = test_stl_desugar
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
% Written:       12-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% setup
x = stl('x', 3);

% cases
p1 = finally(globally(x(1) < 5, interval(1,2)), interval(2,3));
r1 = until(stl(true), ~ until(stl(true), ~(x(1) < 5), interval(1,2)), interval(2,3));

p2 = release(x(1) < 3, x(2) > 7, interval(1,2));
r2 = ~until(~(x(1) < 3), ~(x(2) > 7), interval(1,2));

p3 = implies(x(1) < 4, finally(x(2) > 5, interval(1,2)));
r3 = ~(x(1) < 4) | until(stl(true), x(2) > 5, interval(1,2));

% test
res = isequal(r1, desugar(p1));
res(end+1,1) = res && isequal(r2, desugar(p2));
res(end+1,1) = isequal(r3, desugar(p3));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
