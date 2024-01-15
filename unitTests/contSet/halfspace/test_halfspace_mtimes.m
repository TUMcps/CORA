function res = test_halfspace_mtimes
% test_halfspace_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_halfspace_mtimes
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
% Written:       16-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
h = halfspace.empty(2);
A = [1 0; 0 2];
res(end+1,1) = representsa(A*h,'emptySet');

% 2D halfspace
c = [1;1]; d = 0.5;
h = halfspace(c,d);
A = [2 3; 1 2];
h_mapped = A*h;
% compare to true result
A_inv = [2 -3; -1 2];
c_true = A_inv' * c;
h_true = halfspace(c_true,d);
res(end+1,1) = isequal(h_mapped,h_true);

% combine results
res = all(res);

% test singular matrix multiplication
try 
    A = [1 0; 0 0];
    A * h; % should throw an error
    res = false;
    return
end

% ------------------------------ END OF CODE ------------------------------
