function res = test_specification_isempty
% test_specification_isempty - unit test for isempty
%
% Syntax:
%    res = test_specification_isempty
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
% Written:       02-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init helpers
set1 = zonotope([0;0],[1,-0.7;0.2,1]);
set2 = interval([-1;-2],[2;3]);
set3 = zonotope.empty(2);
time1 = interval(2,4);
time2 = interval(1,3);

spec1 = specification(set1,'unsafeSet',time1);
spec2 = specification(set2,'safeSet',time2);
spec3 = specification(set3,'invariant');
spec12 = [spec1;spec2];
spec13 = [spec1;spec3];

spec_empty = specification();

% check emptiness
assert(isempty(spec_empty));
assert(~isempty(spec1));
assert(~isempty(spec2));
assert(~isempty(spec3));
assert(~isempty(spec12));
assert(~isempty(spec13));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
