function res = test_interval_vertices
% test_interval_vertices - unit test function of vertices
%
% Syntax:  
%    res = test_interval_vertices
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      28-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty set
I_e = interval();
res = isnumeric(vertices(I_e)) && isempty(vertices(I_e));

% create interval
lowerLimits = [-2; -4];
upperLimits = [3; 1];
int = interval(lowerLimits, upperLimits);

% retreive vertices by function
vert = vertices(int);
vert_mat = vert;
vert_mat = sortrows(vert_mat');

% true vertices
vert_true = [-2 3  3 -2;
              1 1 -4 -4];
vert_true = sortrows(vert_true');

% compare results
tol = 1e-9;
res = res && all(all(abs(vert_mat - vert_true) < tol));

%------------- END OF CODE --------------