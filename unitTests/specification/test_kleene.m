function res = test_kleene
% test_kleene - unit test function of kleene
%
% Syntax:
%    res = test_kleene
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
% Written:       10-August-2022
% Last update:   21-February-2024 (FL, add tests for fromBool and array operators)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% conjunction
assert(isequal(kleene.True, kleene.True & kleene.True));
assert(isequal(kleene.False, kleene.True & kleene.False));
assert(isequal(kleene.Unknown, kleene.True & kleene.Unknown));
assert(isequal(kleene.False, kleene.False & kleene.True));
assert(isequal(kleene.False, kleene.False & kleene.False));
assert(isequal(kleene.False, kleene.False & kleene.Unknown));
assert(isequal(kleene.Unknown, kleene.Unknown & kleene.True));
assert(isequal(kleene.False, kleene.Unknown & kleene.False));
assert(isequal(kleene.Unknown, kleene.Unknown & kleene.Unknown));

% disjunction
assert(isequal(kleene.True, kleene.True | kleene.True));
assert(isequal(kleene.True, kleene.True | kleene.False));
assert(isequal(kleene.True, kleene.True | kleene.Unknown));
assert(isequal(kleene.True, kleene.False | kleene.True));
assert(isequal(kleene.False, kleene.False | kleene.False));
assert(isequal(kleene.Unknown, kleene.False | kleene.Unknown));
assert(isequal(kleene.True, kleene.Unknown | kleene.True));
assert(isequal(kleene.Unknown, kleene.Unknown | kleene.False));
assert(isequal(kleene.Unknown, kleene.Unknown | kleene.Unknown));

% negation
assert(isequal(kleene.False, ~ kleene.True));
assert(isequal(kleene.True, ~ kleene.False));
assert(isequal(kleene.Unknown, ~ kleene.Unknown));

% fromBool
assert(isequal(kleene.True, kleene.fromBool(true)));
assert(isequal(kleene.False, kleene.fromBool(false)));

% arrays
arr1 = [kleene.True, kleene.False, kleene.Unknown];
arr2 = [kleene.Unknown, kleene.False, kleene.True];
assert(isequal([kleene.Unknown, kleene.False, kleene.Unknown], arr1 & arr2));
assert(isequal([kleene.True, kleene.False, kleene.True], arr1 | arr2));
assert(isequal([kleene.False, kleene.True, kleene.Unknown], ~arr1));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
