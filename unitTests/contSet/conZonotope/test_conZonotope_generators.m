function res = test_conZonotope_generators
% test_conZonotope_generators - unit test function for reading out the
%    generator matrix of a conZonotope object
%
% Syntax:
%    res = test_conZonotope_generators
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Mark Wetzlinger
% Written:       28-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty conZonotope
n = 2;
cZ_empty = conZonotope.empty(n);
G =  generators(cZ_empty);
assert(isempty(G) && isnumeric(G) && size(G,1) == n);

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
