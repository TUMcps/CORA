function [Zbundle] = quadraticMultiplication_zono(Zbundle,Q)
% quadraticMultiplication_zono - computes \{Q_{ijk}*x_j*x_k|x \in Z\} when
% Q is modeled as a tensor zonotope
%
% Syntax:  
%    [Zquad] = quadraticMultiplication_zono(Zbundle,Q)
%
% Inputs:
%    Zbundle - zonotopeBundle object
%    Q - quadratic coefficients as a cell of matrix zonotopes
%
% Outputs:
%    Zquad - zonotopeBundle object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      26-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%compute quadratic map for each zonotope pair
for i=1:Zbundle.parallelSets
    Zbundle.Z{i} = quadraticMultiplication_zono(Zbundle.Z{i},Q);
end


%------------- END OF CODE --------------