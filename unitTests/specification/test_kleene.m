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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% conjunction
res(end+1,1) = isequal(kleene.True, kleene.True & kleene.True);
res(end+1,1) = isequal(kleene.False, kleene.True & kleene.False);
res(end+1,1) = isequal(kleene.Unknown, kleene.True & kleene.Unknown);
res(end+1,1) = isequal(kleene.False, kleene.False & kleene.True);
res(end+1,1) = isequal(kleene.False, kleene.False & kleene.False);
res(end+1,1) = isequal(kleene.False, kleene.False & kleene.Unknown);
res(end+1,1) = isequal(kleene.Unknown, kleene.Unknown & kleene.True);
res(end+1,1) = isequal(kleene.False, kleene.Unknown & kleene.False);
res(end+1,1) = isequal(kleene.Unknown, kleene.Unknown & kleene.Unknown);

% disjunction
res(end+1,1) = isequal(kleene.True, kleene.True | kleene.True);
res(end+1,1) = isequal(kleene.True, kleene.True | kleene.False);
res(end+1,1) = isequal(kleene.True, kleene.True | kleene.Unknown);
res(end+1,1) = isequal(kleene.True, kleene.False | kleene.True);
res(end+1,1) = isequal(kleene.False, kleene.False | kleene.False);
res(end+1,1) = isequal(kleene.Unknown, kleene.False | kleene.Unknown);
res(end+1,1) = isequal(kleene.True, kleene.Unknown | kleene.True);
res(end+1,1) = isequal(kleene.Unknown, kleene.Unknown | kleene.False);
res(end+1,1) = isequal(kleene.Unknown, kleene.Unknown | kleene.Unknown);

% negation
res(end+1,1) = isequal(kleene.False, ~ kleene.True);
res(end+1,1) = isequal(kleene.True, ~ kleene.False);
res(end+1,1) = isequal(kleene.Unknown, ~ kleene.Unknown);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
