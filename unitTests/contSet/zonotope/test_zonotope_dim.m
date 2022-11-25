function res = test_zonotope_dim
% test_zonotope_dim - unit test function of dim
%
% Syntax:  
%    res = test_zonotope_dim
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% check empty zonotope
Z = zonotope();
res(1) = dim(Z) == 0;


% instantiate zonotope
c = [-2; 1];
G = [2 4 5 3 3; 0 3 5 2 3];
Z = zonotope(c,G);

res(2) = dim(Z) == 2;

% no generator matrix
Z = zonotope(c);
res(3) = dim(Z) == 2;

% combine results
res = all(res);

if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------