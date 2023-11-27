function res = test_removeDuplicates()
% test_removeDuplicates - unit test function for the removal of duplicate
%    vertices
%
% Syntax:
%    res = test_removeDuplicates()
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
% Written:       21-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% too few vertices
V = [1 0]';
V_ = removeDuplicates(V);
res(end+1,1) = compareMatrices(V,V_);

% init vertices
V = [-1 0; 0 0; 1 0; 1 0]';
V_ = removeDuplicates(V);
res(end+1,1) = compareMatrices(V_,[-1 0; 0 0; 1 0]');

% init vertices
V = [-1 0 1; 0 0 1; 1 0 1; 0 0 1; 1 0 1; 2 0 0; 0 0 1]';
V_ = removeDuplicates(V);
res(end+1,1) = compareMatrices(V_,[-1 0 1; 0 0 1; 1 0 1; 2 0 0]');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
