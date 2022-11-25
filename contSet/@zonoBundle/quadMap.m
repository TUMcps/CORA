function [Z] = quadMap(Z,Q)
% quadMap - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:  
%    [Zquad] = quadMap(Z1,Q)
%
% Inputs:
%    Z - zonoBundle object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    Z - zonoBundle object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap

% Author:       Niklas Kochdumper
% Written:      13-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    for i=1:Z.parallelSets
        Z.Z{i} = quadMap(Z.Z{i},Q);
    end

%------------- END OF CODE --------------