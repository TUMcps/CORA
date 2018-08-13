function [Zquad] = quadraticMultiplication(Z,Q)
% quadraticMultiplication - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:  
%    [Zquad] = quadraticMultiplication(Z1,Q)
%
% Inputs:
%    Z - zonotope object
%    Q - quadratic coefficients as a cell of matrices
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
% Written:      07-December-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%get matrix of zonotope
Zmat = get(Z,'Z');
dim = length(Q);
gens = length(Zmat(1,:)) - 1;

%for each dimension, compute generator elements
for i = 1:dim
    
    %pure quadratic evaluation
    quadMat = Zmat'*Q{i}*Zmat;
    
    %center
    c(i,1) = quadMat(1,1) + 0.5*sum(diag(quadMat(2:end,2:end)));
    
    %generators from diagonal elements
    ind = 1:gens;
    G(i, ind) = 0.5*diag(quadMat(ind+1,ind+1));
    
    %generators from other elements
    counter = 0;
    for j = 0:gens
        kInd = (j+1):gens;
        G(i, gens + counter + kInd - j) = quadMat(j+1, kInd+1) + quadMat(kInd+1, j+1)';
        counter = counter + length(kInd);
    end
    
end

%generate new zonotope
Zquad = zonotope([c, G]);


%------------- END OF CODE --------------