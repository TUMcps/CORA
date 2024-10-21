function res = test_finiteKleeneSignalBuilder
% test_finiteKleeneSignalBuilder - unit test of finiteKleeneSignalBuilder class
%
% Syntax:
%    res = test_finiteKleeneSignalBuilder
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
% Written:       23-August-2022
% Last update:   14-February-2024 (FL, rename kleeneSignalBuilder to finiteKleeneSignalBuilder)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% builders
b1 = finiteKleeneSignalBuilder(10);
b2 = finiteKleeneSignalBuilder(10);
b3 = finiteKleeneSignalBuilder(10);

% test
b1.setPossiblyTrue(interval(0, 5));
b1.setPossiblyFalse(interval(3, 7));
b1.setPossiblyTrue(interval(7, 9));
b1.setPossiblyTrue(interval(8, 10));
b1.setPossiblyFalse(interval(9, 10));

s1 = b1.kleeneSignal();

assert(isequal([3 5 7 9 10], s1.time));
assert(isequal([tt uu ff tt uu], s1.value));

% test
b2.setPossiblyTrue(interval(0, 10));
b2.setPossiblyFalse(interval(0, 10));

s2 = b2.kleeneSignal();

assert(isequal(10, s2.time));
assert(isequal(uu, s2.value));

% test
b3.setPossiblyFalse(interval(0, 2));
b3.setPossiblyTrue(interval(2, 4));
b3.setPossiblyTrue(interval(4, 6));
b3.setPossiblyFalse(interval(6, 8));
b3.setPossiblyTrue(interval(8, 10));

s3 = b3.kleeneSignal();

assert(isequal([2 6 8 10], s3.time));
assert(isequal([ff tt ff tt], s3.value));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
