function res = test_nonzeroFilter()
% test_nonzeroFilter - unit test function for deletion of 0-length
%    generators
%
% Syntax:
%    res = test_nonzeroFilter()
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
% Written:       24-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% matrix without zeros
G = [1 2 3; 4 5 6];
G_filtered = nonzeroFilter(G);
res(end+1,1) = compareMatrices(G_filtered,G);

% matrix with zeros
G = [1 2 0 3 0; 4 5 0 6 0];
G_filtered = nonzeroFilter(G);
G_true = [1 2 3; 4 5 6];
res(end+1,1) = compareMatrices(G_filtered,G_true);

% matrix with almost zeros, zero tolerance
G = [1 2 eps 3; 4 5 eps 6];
G_filtered = nonzeroFilter(G);
G_true = [1 2 3; 4 5 6];
res(end+1,1) = ~compareMatrices(G_filtered,G_true);

% matrix with almost zeros, enough tolerance
G = [1 2 eps 3; 4 5 eps 6];
tol = 2*eps;
G_filtered = nonzeroFilter(G,tol);
G_true = [1 2 3; 4 5 6];
res(end+1,1) = compareMatrices(G_filtered,G_true);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
