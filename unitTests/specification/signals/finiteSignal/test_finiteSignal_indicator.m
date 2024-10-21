function res = test_finiteSignal_indicator
% test_finiteSignal_indicator - unit test of signal indicator function
%
% Syntax:
%    res = test_finiteSignal_indicator
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
% Written:       16-August-2022
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% signals
sig1 = finiteSignal.indicator(3.0, interval(2.3, 2.7), tt);
sig2 = finiteSignal.indicator(3.0, interval(2.3, 2.7), uu);
sig3 = finiteSignal.indicator(3.0, interval(2.3, 2.7), ff);

sig4 = finiteSignal.indicator(3.0, interval(0.0, 2.7), tt);
sig5 = finiteSignal.indicator(3.0, interval(2.3, 3.0), tt);
sig6 = finiteSignal.indicator(3.0, interval(0.0, 3.0), tt);

sig7 = finiteSignal.indicator(3.0, interval(2.7, 4.3), tt);
sig8 = finiteSignal.indicator(3.0, interval(0, 4.3), tt);
sig9 = finiteSignal.indicator(3.0, interval(-1.3, 1.4), tt);

sig10 = finiteSignal.indicator(3.0, interval(4.2, 5.7), tt);
sig11 = finiteSignal.indicator(3.0, interval(-2.7, -2.4), tt);
sig12 = finiteSignal.indicator(3.0, interval(-2.3, 3.7), tt);

% test
assert(isequal([2.3 2.7 3.0], sig1.time));
assert(isequal([ff tt ff], sig1.value));

assert(isequal([2.3, 2.7, 3.0], sig2.time));
assert(isequal([ff uu ff], sig2.value));

assert(isequal(3.0, sig3.time));
assert(isequal(ff, sig3.value));

assert(isequal([2.7 3.0], sig4.time));
assert(isequal([tt ff], sig4.value));

assert(isequal([2.3 3.0], sig5.time));
assert(isequal([ff tt], sig5.value));

assert(isequal(3.0, sig6.time));
assert(isequal(tt, sig6.value));

assert(isequal([2.7 3.0], sig7.time));
assert(isequal([ff tt], sig7.value));

assert(isequal(3.0, sig8.time));
assert(isequal(tt, sig8.value));

assert(isequal([1.4 3.0], sig9.time));
assert(isequal([tt ff], sig9.value));

assert(isequal(3.0, sig10.time));
assert(isequal(ff, sig10.value));

assert(isequal(3.0, sig11.time));
assert(isequal(ff, sig11.value));

assert(isequal(3.0, sig12.time));
assert(isequal(tt, sig12.value));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
