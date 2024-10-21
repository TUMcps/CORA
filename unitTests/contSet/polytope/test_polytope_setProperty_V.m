function res = test_polytope_setProperty_V
% test_polytope_setProperty_V - unit test function to check whether the
%    internally-used set property 'V' is changed correctly following
%    different set operations on a polytope
%
% Syntax:
%    res = test_polytope_setProperty_V
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
% Written:       14-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% --- mtimes --------------------------------------------------------------

% 2D, origin
V = [0;0];
P = polytope(V);
A = [2 0; -1 1];
P_ = A * P;
assert(P_.isVRep.val && compareMatrices(P_.V,V));

% 2D, bounded
P = polytope([1 0; -1 1; -1 -1], [1; 1; 1]);
M = [1 -1; 2 1];
V = vertices(P);
% ... vertices now stored in P
P_ = M*P;
vertices(P_);
% ... vertices now stored in P_
assert(P.isVRep.val && compareMatrices(P.V,V));
assert(P_.isVRep.val && compareMatrices(P_.V,M*V));

% --- minus ---------------------------------------------------------------

% 1D, bounded
V = [-0.5, 1];
P = polytope(V);
P_ = P - 2;
V_ = V - 2;
assert(P_.isVRep.val && compareMatrices(P_.V,V_));

% 2D, bounded
V = [1 0; -1 1; 0 -1]';
P = polytope(V);
P_ = P - [1;-1];
V_ = V - [1;-1];
assert(P_.isVRep.val && compareMatrices(P_.V,V_));

% --- plus ----------------------------------------------------------------

% 1D, bounded
V = [-0.5, 1];
P = polytope(V);
P_ = P + 2;
V_ = V + 2;
assert(P_.isVRep.val && compareMatrices(P_.V,V_));

% 2D, bounded
V = [1 0; -1 1; 0 -1]';
P = polytope(V);
P_ = P + [1;-1];
V_ = V + [1;-1];
assert(P_.isVRep.val && compareMatrices(P_.V,V_));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
