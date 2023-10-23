function res = test_zonoBundle_contains
% test_zonoBundle_contains - unit test function of contains
%
% Syntax:
%    res = test_zonoBundle_contains
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

% Authors:       Mark Wetzlinger
% Written:       16-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate zonotope bundle
Z1 = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
Z2 = Z1 + [1;0;0];
zB = zonoBundle({Z1,Z2});

% points inside/outside
p_in = [0;0;0];
p_out = [4;0.5;1];

% check points
res = contains(zB,p_in);
res(end+1) = contains(zB,p_in,'exact:zonotope');
res(end+1) = contains(zB,p_in,'exact:polytope');
res(end+1) = ~contains(zB,p_out);
res(end+1) = ~contains(zB,p_out,'exact:zonotope');
res(end+1) = ~contains(zB,p_out,'exact:polytope');

% interval inside
I = interval(-0.02*ones(3,1),0.02*ones(3,1));
res(end+1) = contains(zB,I);
res(end+1) = contains(zB,I,'exact:polytope');

% zonotope inside
Z = zonotope(zeros(3,1),0.02*[1 2 0; -2 1 1; 0 -1 1]);
res(end+1) = contains(zB,Z);

% capsule inside
C = capsule(zeros(3,1),0.02*ones(3,1),0.02);
res(end+1) = contains(zB,C);

% ellipsoid inside
E = ellipsoid(0.001*[3.2 1.1 -0.3; 1.1 11.4 3.9; -0.3 3.9 2.6],zeros(3,1));
res(end+1) = contains(zB,E);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
