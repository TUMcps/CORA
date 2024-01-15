function res = test_halfspace_plus
% test_halfspace_plus - unit test function of plus
%
% Syntax:
%    res = test_halfspace_plus
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
% Written:       16-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
h = halfspace.empty(2);
v = [2; 1];
res(end+1,1) = representsa(h + v,'emptySet');

% combine results
res = all(res);


% dimension mismatch
h = halfspace(randn(2,1),1);
v = rand(3,1);
try
    h + v; % should throw error here
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
