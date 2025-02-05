function res = test_finiteSignal_cutoff
% test_finiteSignal_cutoff - unit test of signal cutoff function
%
% Syntax:
%    res = test_finiteSignal_cutoff
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
% Written:       19-August-2022
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% signals
sig1 = finiteSignal([1.2 2.3 3.0], [true false true]);
sig2 = finiteSignal(2.7, false);

% test
c1 = cutoff(sig1, 2.7);

assert(isequal([1.2 2.3 2.7], c1.time));
assert(isequal([true false true], c1.value));

c2 = cutoff(sig1, 3.0);

assert(isequal([1.2 2.3 3.0], c2.time));
assert(isequal([true false true], c2.value));

c3 = cutoff(sig1, 2.0);

assert(isequal([1.2 2.0], c3.time));
assert(isequal([true false], c3.value));

c4 = cutoff(sig1, 2.3);

assert(isequal([1.2 2.3], c4.time));
assert(isequal([true false], c4.value));

% more tests
c5 = cutoff(sig1, 0.1);

assert(isequal(0.1, c5.time));
assert(isequal(true, c5.value));

c6 = cutoff(sig1, 0);

assert(isequal(0, c6.time));
assert(isequal(true, c6.value));

c7 = cutoff(sig2, 2.7);

assert(isequal(2.7, c7.time));
assert(isequal(false, c7.value));

c8 = cutoff(sig2, 1.2);

assert(isequal(1.2, c8.time));
assert(isequal(false, c8.value));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
