function res = test_checkValueAttributes
% test_checkValueAttributes - unit test function for checkValueAttributes
%
% Syntax:
%    res = test_checkValueAttributes()
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
% See also: checkValueAttributes

% Authors:       Tobias Ladner
% Written:       03-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% numeric -----------------------------------------------------------------

% check class
assert(checkValueAttributes(1,'numeric',{}))

% check attribute properties ---

% check single
assert(checkValueAttributes(1,'numeric',{'positive'}))
assert(checkValueAttributes(1,'numeric',{'nonnegative'}))
assert(checkValueAttributes(1,'numeric',{'scalar'}))
assert(checkValueAttributes(1:2,'numeric',{@(dims) all(dims == [1,2])}))

% check multiple
assert(checkValueAttributes(1,'numeric',{'integer','positive'}))

% check unknown
assertThrowsAs(@checkValueAttributes,'MATLAB:UndefinedFunction',1,'numeric',{'unknown'})

% char --------------------------------------------------------------------

% check class
assert(checkValueAttributes('test','char',{}))

% check single
assert(checkValueAttributes('mystring','char',{'nonempty'}))

% contSet -----------------------------------------------------------------

Z = zonotope([1;2],[1 0; 0 1]);
assert(checkValueAttributes(Z,'contSet',{@(Z) representsa(Z,'interval')}))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
