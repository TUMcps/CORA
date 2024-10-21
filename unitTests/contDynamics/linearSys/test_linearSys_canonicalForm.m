function res = test_linearSys_canonicalForm
% test_linearSys_canonicalForm - unit test for the computation of the
%    canonical form of linear systems
%
% Syntax:
%    res = test_linearSys_canonicalForm
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-10;

% init linear system
A = [-1 -4; 4 -1];
B = [1; 0];
c = [4; -2];
C = [1 0; 0 -1; 3 2];
D = [1; -1; 4];
k = [4; 2; 0];
E = eye(2);
F = eye(3);
sys = linearSys(A,B,c,C,D,k,E,F);

% init sets
U = zonotope(1,2);
uVec = [1 2 1 0 -1];
W = zonotope([1;0],[1 -2 1; 0 3 1]);
V = zonotope([1;-1;2], [1 0 0; 0 1 2; 0 0 -1]);
vVec = [1 0 0; 1 1 0; 1 2 0; 1 2 1; 1 2 2; 1 2 1]';

% compute canonical form
% Ax + Bu + c + w  ->  Ax + u', where u' \in U_ + u_
% Cx + Du + k + v  ->  Cx + v', where v' \in V_ + v_
[sys,U_,u_,V_,v_] = canonicalForm(sys,U,uVec,W,V,vVec);

steps = size(uVec,2);
U_true = B*U + E*W + c;
u_true = B*uVec;
V_true = D*U + F*V + k;
v_true = D*[uVec, 0] + vVec;
Uu_true = arrayfun(@(i) U_true+u_true(:,i), 1:steps, 'UniformOutput', false);
Vv_true = arrayfun(@(i) V_true+v_true(:,i), 1:steps, 'UniformOutput', false);

% it suffices to ensure that the sets (U_+u_) and (V_+v_) are correct
for i=1:steps
    assert(isequal(U_+u_(:,i),Uu_true{i},tol));
    assert(isequal(V_+v_(:,i),Vv_true{i},tol));
end


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
