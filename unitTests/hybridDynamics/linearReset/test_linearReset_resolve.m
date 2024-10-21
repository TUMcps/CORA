function res = test_linearReset_resolve
% test_linearReset_resolve - test function for the resolution of local
%    inputs to states of other components (via outputs) and global inputs
%
% Syntax:
%    res = test_linearReset_resolve
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
% See also: test_nonlinearReset_resolve

% Authors:       Mark Wetzlinger
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = eps;

% init synchronized linear reset function (original values for A and c are
% only constant offset to resulting values, so we can keep them at zero)
A = zeros(5);
B1 = [1 -1; 2 -2; 3 -3]; B2 = [4; 5];
B = blkdiag(B1, B2);
c = zeros(5,1);
linReset = linearReset(A,B,c);

% binds and flows (state equation does not matter here)
stateBinds = {[1,2,3],[4,5]};
inputBinds = {[2 1; 0 1],[1 1]};
sys1 = linearSys(eye(3),zeros(3,2),zeros(3,1),[1,2,3],[0,100],5);
sys2 = linearSys(eye(2),[0;0],[0;0],[-1,-2],0,-10);
flowList = {sys1;sys2};

% resolve input binds
linReset_ = resolve(linReset,flowList,stateBinds,inputBinds);
% A: affected by resolution of u11 = y21 and u21 = y11
A_resolved = [zeros(3,3), B1(:,1) * sys2.C; B2(:,1) * sys1.C, zeros(2,2)];
assert(compareMatrices(linReset_.A,A_resolved,tol,"equal",true));
% B: only part with global inputs remains
B_resolved = [B1(:,2); B2*sys1.D(2)];
assert(compareMatrices(linReset_.B,B_resolved,tol,"equal",true));
% c: affected by resolution of u11 = y21 and u21 = y11
c_resolved = [B1(:,1) * sys2.k; B2(:,1) * sys1.k];
assert(compareMatrices(linReset_.c,c_resolved,tol,"equal",true));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
