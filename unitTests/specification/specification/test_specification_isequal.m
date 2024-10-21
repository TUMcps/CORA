function res = test_specification_isequal
% test_specification_isequal - unit test for isequal
%
% Syntax:
%    res = test_specification_isequal
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
assert(isequal(spec1,spec1));
assert(~isequal(spec1,spec2));

% with time
spec1 = specification(I1,'safeSet',time1);
spec2 = specification(I1,'safeSet',time2);
spec3 = specification(I1,'safeSet');
assert(isequal(spec1,spec1));
assert(~isequal(spec1,spec2));
assert(~isequal(spec1,spec3));

% with location
spec1 = specification(I1,'unsafeSet',loc1);
spec2 = specification(I1,'unsafeSet',loc2);
spec3 = specification(I1,'unsafeSet');
assert(isequal(spec1,spec1));
assert(~isequal(spec1,spec2));
assert(~isequal(spec1,spec3));

% with tolerance
spec1 = specification(I1);
spec2 = specification(I1+eps);
assert(isequal(spec1,spec2));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
