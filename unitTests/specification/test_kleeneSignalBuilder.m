function res = test_kleeneSignalBuilder
% test_kleeneSignalBuilder - unit test of kleeneSignalBuilder class
%
% Syntax:
%    res = test_kleeneSignalBuilder
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% shortcodes
tt = kleene.True;
ff = kleene.False;
uu = kleene.Unknown;

% builders
b1 = kleeneSignalBuilder(10);
b2 = kleeneSignalBuilder(10);
b3 = kleeneSignalBuilder(10);

% test
b1.setPossiblyTrue(interval(0, 5));
b1.setPossiblyFalse(interval(3, 7));
b1.setPossiblyTrue(interval(7, 9));
b1.setPossiblyTrue(interval(8, 10));
b1.setPossiblyFalse(interval(9, 10));

s1 = b1.kleeneSignal();

res(end+1,1) = isequal([3 5 7 9 10], s1.time);
res(end+1,1) = isequal([tt uu ff tt uu], s1.value);

% test
b2.setPossiblyTrue(interval(0, 10));
b2.setPossiblyFalse(interval(0, 10));

s2 = b2.kleeneSignal();

res(end+1,1) = isequal(10, s2.time);
res(end+1,1) = isequal(uu, s2.value);

% test
b3.setPossiblyFalse(interval(0, 2));
b3.setPossiblyTrue(interval(2, 4));
b3.setPossiblyTrue(interval(4, 6));
b3.setPossiblyFalse(interval(6, 8));
b3.setPossiblyTrue(interval(8, 10));

s3 = b3.kleeneSignal();

res(end+1,1) = isequal([2 6 8 10], s3.time);
res(end+1,1) = isequal([ff tt ff tt], s3.value);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
