function res = test_conZonotope_interval
% test_conZonotope_interval - unit test function for the caclulation of
%                             a bounding box of a constrained zonotope object
%
% Syntax:
%    res = test_conZonotope_interval
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

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       22-May-2018
% Last update:   23-February-2024 (MW, add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% % Figure 1 in [1]
% Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
% A = [1 1 1]; b = 1;
% cZ = conZonotope(Z,A,b);
% I = interval(cZ);
% I_ = interval([-2.5;-1.5],[3.5;2.5]);
% res(end+1,1) = isequal(I,I_);
% 
% 
% % Figure 2 in [1]
% Z = [0 1 0 1;0 1 2 -1];
% A = [-2 1 -1]; b = 2;
% cZ = conZonotope(Z,A,b);
% I = interval(cZ);
% I_ = interval([-2;-2],[0;3]);
% res(end+1,1) = isequal(I,I_);


% empty case: constraint beta = -2 infeasible within |beta| <= 1
c = [0; 0];
G = [1; -1];
A = 1; b = -2;
cZ = conZonotope(c,G,A,b);
I = interval(cZ);
assert(representsa_(I,'emptySet',eps));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
