function res = test_specification_ne
% test_specification_ne - unit test for '~=' operator
%
% Syntax:
%    res = test_specification_ne
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init helpers
I1 = interval([-2;1],[3;2]);
I2 = interval([-2;4],[3;5]);
time1 = interval(0,2);
time2 = interval(1,2);
loc1 = [1 2];
loc2 = [1 3];

% with types
spec1 = specification(I1,'safeSet');
spec2 = specification(I2,'unsafeSet');
assert(~(spec1 ~= spec1));
assert(ne(spec1,spec2));

% with time
spec1 = specification(I1,'safeSet',time1);
spec2 = specification(I1,'safeSet',time2);
spec3 = specification(I1,'safeSet');
assert(~(spec1 ~= spec1));
assert(ne(spec1,spec2));
assert(spec1 ~= spec3);

% with location
spec1 = specification(I1,'unsafeSet',loc1);
spec2 = specification(I1,'unsafeSet',loc2);
spec3 = specification(I1,'unsafeSet');
assert(~(spec1 ~= spec1));
assert(ne(spec1,spec2));
assert(spec1 ~= spec3);

% with tolerance
spec1 = specification(I1);
spec2 = specification(I1+eps);
assert(~(spec1 ~= spec2));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
