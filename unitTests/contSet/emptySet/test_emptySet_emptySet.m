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

res = true(0);

% dimension is zero
n = 0;
O = emptySet(n);
res(end+1,1) = O.dimension == n;

% dimension greater or equal to one
n = 2;
O = emptySet(n);
res(end+1,1) = O.dimension == n;

% combine results
res = all(res);


% too many input arguments
try
    O = emptySet(n,n);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
