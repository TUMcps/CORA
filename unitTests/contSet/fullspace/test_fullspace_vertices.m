function res = test_fullspace_vertices
% test_fullspace_vertices - unit test function of vertices
%
% Syntax:  
%    res = test_fullspace_vertices
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
% Written:      25-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init fullspace
n = 2;
fs = fullspace(n);

% compute vertices
V = vertices(fs);

% check result
res = all(contains(fs,V));

%------------- END OF CODE --------------