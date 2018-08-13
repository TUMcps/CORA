function [Zquad] = quadraticMultiplication_zono(Z,Q)
% quadraticMultiplication_zono - computes \{Q_{ijk}*x_j*x_k|x \in Z\} when
% Q is modeled as a tensor zonotope
%
% Syntax:  
%    [Zquad] = quadraticMultiplication_zono(Z,Q)
%
% Inputs:
%    Z - zonotope object
%    Q - quadratic coefficients as a cell of matrix zonotopes
%
% Outputs:
%    Zquad - zonotope object
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
% Written:      26-September-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%compute quadratic map of center tensor
parfor i = 1:length(Q)
    C{i} = Q{i}.center;
end
Z_c = quadraticMultiplication(Z,C);


%compute quadratic map of generator tensors
parfor iGen = 1:Q{1}.gens
    for i=1:length(Q)
        G{iGen}{i} = Q{i}.generator{iGen};
    end
end
parfor iGen = 1:Q{1}.gens  
    Z_g{iGen} = quadraticMultiplication(Z,G{iGen}(:));
end

%combine all zonotopes to a single zonotope
%center matrix
Zmat_center = get(Z_c,'Z');

%generator matrices
Zmat_gen = [];
newGens = length(Zmat_center(1,:));
for iGen = 1:Q{1}.gens
   Zmat_gen(:,(end+1):(end+newGens)) = get(Z_g{iGen},'Z');
end

%generate zonotope
Zquad = zonotope([Zmat_center, Zmat_gen]);

%------------- END OF CODE --------------