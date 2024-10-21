function res = test_polytope_origin
% test_polytope_origin - unit test function of origin
%
% Syntax:
%    res = test_polytope_origin
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
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D
P = polytope.origin(1);
P_true = polytope(0);
assert(isequal(P,P_true));
assert(contains(P,0));
assert(~isempty(P.emptySet) && ~P.emptySet.val);
assert(~isempty(P.fullDim) && ~P.fullDim.val);
assert(~isempty(P.bounded) && P.bounded.val);
assert(~isempty(P.minHRep) && P.minHRep.val);
assert(~isempty(P.minVRep) && P.minVRep.val);

% 2D
P = polytope.origin(2);
P_true = polytope(zeros(2,1));
assert(isequal(P,P_true));
assert(contains(P,zeros(2,1)));
assert(~isempty(P.emptySet) && ~P.emptySet.val);
assert(~isempty(P.fullDim) && ~P.fullDim.val);
assert(~isempty(P.bounded) && P.bounded.val);
assert(~isempty(P.minHRep) && P.minHRep.val);
assert(~isempty(P.minVRep) && P.minVRep.val);

% wrong calls
assertThrowsAs(@polytope.origin,'CORA:wrongValue',0);
assertThrowsAs(@polytope.origin,'CORA:wrongValue',-1);
assertThrowsAs(@polytope.origin,'CORA:wrongValue',0.5);
assertThrowsAs(@polytope.origin,'CORA:wrongValue',[1,2]);
assertThrowsAs(@polytope.origin,'CORA:wrongValue','text');

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
