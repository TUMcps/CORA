function res = test_stl_masks
% test_stl_masks - unit test of stl masks method
%
% Syntax:
%    res = test_stl_masks
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
% Written:       23-January-2023
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

x = stl('x', atomicProposition(interval(0,1)));
y = stl('y', atomicProposition(interval(0,1)));
z = stl('z', atomicProposition(interval(0,1)));

phi1 = finally(x, interval(1.2, 1.4));
phi2 = until(x, y & z, interval(2.3, 4));
phi3 = release(phi1, ~z, interval(0, 1.7));
phi4 = release(phi2, ~z, interval(0, 1.7));

m1 = masks(phi1);
m2 = masks(phi2);
m3 = masks(phi3);
m4 = masks(phi4);

% test masks function
res = true;
assert(isequal(finiteSignal([1.2 1.4], [false true]), m1("x")));

assert(isequal(finiteSignal(4, true), m2("x")));
assert(isequal(finiteSignal([2.3 4], [false true]), m2("y")));
assert(isequal(finiteSignal([2.3 4], [false true]), m2("z")));

assert(isequal(finiteSignal([1.2 1.4+1.7], [false true]), m3("x")));
assert(isequal(finiteSignal([1.7 1.4+1.7], [true false]), m3("z")));

assert(isequal(finiteSignal(5.7, true), m4("x")));
assert(isequal(finiteSignal([2.3 5.7], [false true]), m4("y")));
assert(isequal(finiteSignal([1.7 2.3 5.7], [true false true]), m4("z")));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
