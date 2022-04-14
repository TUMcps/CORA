function Zquad = quadMap_parallel(Z,Q)
% quadMap - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:  
%    Zquad = quadMap(Z1,Q)
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
Zmat = Z.Z;
dimQ = length(Q);
gens = length(Zmat(1,:)) - 1;


%for each dimension, compute generator elements
parfor i = 1:dimQ
    
    
    %pure quadratic evaluation
    quadMat = Zmat'*Q{i}*Zmat;
    
    %center
    ind = 1:gens;
    c(i,1) = quadMat(1,1) + 0.5*sum(diag(quadMat(2:end,2:end)));
    
    %generators with center
    G{i}(ind) = quadMat(1,ind+1) + quadMat(ind+1,1)';
    
    %generators from diagonal elements
    ind = 1:gens;
    G{i}(gens + ind) = 0.5*diag(quadMat(ind+1,ind+1));
    
    %generators from other elements
    counter = 0;
    for j = 1:gens
        kInd = (j+1):gens;
        G{i}(2*gens + counter + kInd - j) = quadMat(j+1, kInd+1) + quadMat(kInd+1, j+1)';
        counter = counter + length(kInd);
    end
end

for i = 1:dimQ
    Gmat(i,:) = G{i}; 
end

%generate new zonotope
Zquad = zonotope([c, Gmat]);


%------------- END OF CODE --------------