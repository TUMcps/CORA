function res = test_conZonotope_contains
% test_conZonotope_contains - unit test function for containment check
%
% Syntax:
%    res = test_conZonotope_contains
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty set
% cZ = conZonotope.empty(2);

% 2D constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1]; b = 1;
cZ = conZonotope(Z,A,b);

% points that are known to be contained
p_in = [0 1 0  0 1 -1 -2; 
        0 0 2 -1 1  0 -1];
% points that are known to be outside
p_out = [3  0 1 -2  3;
         1 -2 2  0 -1];

res(end+1,1) = all(contains(cZ,p_in));
res(end+1,1) = ~any(contains(cZ,p_out));

% visualization
% figure; hold on;
% plot(cZ);
% plot(p_in(1,:),p_in(2,:),'.g');
% plot(p_out(1,:),p_out(2,:),'.r');


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
