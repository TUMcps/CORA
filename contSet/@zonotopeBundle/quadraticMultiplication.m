function [Z] = quadraticMultiplication(Z,Q)
% quadraticMultiplication - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:  
%    [Zquad] = quadraticMultiplication(Z1,Q)
%
% Inputs:
%    Z - zonotopeBundle object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    Z - zonotopeBundle object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      13-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    for i=1:Z.parallelSets
        Z.Z{i} = quadraticMultiplication(Z.Z{i},Q);
    end

%------------- END OF CODE --------------