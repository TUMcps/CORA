function res = test_emptySet_emptySet
% test_emptySet_emptySet - unit test function of constructor
%
% Syntax:
%    res = test_emptySet_emptySet
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension is zero
n = 0;
O = emptySet(n);
assert(O.dimension == n);

% dimension greater or equal to one
n = 2;
O = emptySet(n);
assert(O.dimension == n);

% combine results
res = true;


% too many input arguments
assertThrowsAs(@emptySet,'CORA:numInputArgsConstructor',n,n);

% ------------------------------ END OF CODE ------------------------------
