function res = test_zonotope_isFullDim
% test_zonotope_isFullDim - unit test function of isFullDim
%
% Syntax:  
%    res = test_zonotope_isFullDim
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

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
Z = zonotope();
res(1) = ~isFullDim(Z);

% full-rank matrix
c = [0; 0];
G = [1 2; 2 1];
Z = zonotope(c,G);
res(2) = isFullDim(Z);

% degenerate zonotope
G = [1 0; 0 0];
Z = zonotope(c,G);
res(3) = ~isFullDim(Z);

% combine results
res = all(res);

%------------- END OF CODE --------------