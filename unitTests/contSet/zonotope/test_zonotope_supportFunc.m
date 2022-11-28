function res = test_zonotope_supportFunc
% test_zonotope_supportFunc - unit test function of support function
%
% Syntax:  
%    res = test_zonotope_supportFunc
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

% Author:       Mark Wetzlinger, Victor Gassmann
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;
% empty set
Z_e = polyZonotope();
if supportFunc(Z_e,[1;1],'upper') == -Inf && supportFunc(Z_e,[1;1],'lower') == Inf
    res = true;
end
TOL = 1e-8;

% instantiate zonotope
c = [2; 1];
G = [4 2 -2; 1 3 7];
Z = zonotope(c,G);

%check if [~,x]=minnorm(Z)^2==supportFunc(Z,x)
%check if [~,x]=norm(Z,2,'exact')^2==supportFunc(Z,x)
[val_min,x_min] = minnorm(Z);
[val_max,x_max] = norm(Z,2,'exact');
if abs(supportFunc(Z,x_min)-val_min^2)>TOL || ...
        abs(supportFunc(Z,x_max)-val_max^2)>TOL
    res = true;
end

%------------- END OF CODE --------------