function Z = cubMap(Z,T)
% cubMap - computes an enclosure of the set corresponding to the cubic 
%          multiplication of a zonotope bundle with a third-order tensor
%
% Syntax:  
%    Z = cubMap(Z,Q)
%
% Inputs:
%    Z - zonoBundle object
%    T - third-order tensor
%
% Outputs:
%    Z - zonoBundle object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cubMap

% Author:       Mark Wetzlinger
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

for i=1:Z.parallelSets
    Z.Z{i} = cubMap(Z.Z{i},T);
end

%------------- END OF CODE --------------