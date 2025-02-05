function res = test_finiteSignal_set
% test_finiteSignal_set - unit test of signal set function
%
% Syntax:
%    res = test_finiteSignal_set
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

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
bsig = finiteSignal([1.2 2.3 3.0], [true false true]);
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;
ksig = finiteSignal([1.2 2.3 3.0], [uu ff tt]);

% test
res = true;

sig = bsig.set(stlInterval(1,2), true);
assert(isequal([2.0 2.3 3.0], sig.time));
assert(isequal([true false true], sig.value));

sig = bsig.set(stlInterval(1,2), false);
assert(isequal([1.0 2.3 3.0], sig.time));
assert(isequal([true false true], sig.value));

sig = bsig.set(stlInterval(1,3), false);
assert(isequal([1.0 3.0], sig.time));
assert(isequal([true false], sig.value));

sig = bsig.set(stlInterval(1,3), true);
assert(isequal(3.0, sig.time));
assert(isequal(true, sig.value));

sig = bsig.set(stlInterval(0,1.5), false);
assert(isequal([2.3 3.0], sig.time));
assert(isequal([false true], sig.value));

% more tests
sig = bsig.set(stlInterval(1.2,2.3), true);
assert(isequal(3.0, sig.time));
assert(isequal(true, sig.value));

sig = bsig.set(stlInterval(1.2,2.3), false);
assert(isequal([1.2 2.3 3.0], sig.time));
assert(isequal([true false true], sig.value));

sig = bsig.set(stlInterval(1.2,1.8), true);
assert(isequal([1.8 2.3 3.0], sig.time));
assert(isequal([true false true], sig.value));

sig = ksig.set(stlInterval(1.2,1.8), tt);
assert(isequal([1.2 1.8 2.3 3.0], sig.time));
assert(isequal([uu tt ff tt], sig.value));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
