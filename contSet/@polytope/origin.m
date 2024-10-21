function P = origin(n)
% origin - instantiates a polytope that contains only the origin
%
% Syntax:
%    P = polytope.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    P - polytope representing the origin
%
% Example: 
%    P = polytope.origin(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{n,'att','numeric',{'scalar','positive','integer'}}});

% init halfspace representation (simplex with zero offset)
P = polytope([eye(n); -ones(1,n)],zeros(n+1,1));
% init vertex representation
P.V_.val = zeros(n,1);
P.isVRep.val = true;

% set properties
P.emptySet.val = false;
P.fullDim.val = false;
P.bounded.val = true;
P.minHRep.val = true;
P.minVRep.val = true;

% ------------------------------ END OF CODE ------------------------------
