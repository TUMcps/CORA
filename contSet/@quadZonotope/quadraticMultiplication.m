function [qZquad] = quadraticMultiplication(qZ,Q)
% quadraticMultiplication - computes \{Q_{ijk}*x_j*x_k|x \in qZ\}
%
% Syntax:  
%    [qZquad] = quadraticMultiplication(qZ,Q)
%
% Inputs:
%    qZ - quadZonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    qZquad - quadZonotope object
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
% Written:      05-September-2012
% Last update:  12-August-2017
% Last revision:---

%------------- BEGIN CODE --------------

%get number of dependent  generators
depGens = length(qZ.G(1,:));

%get matrix of overapproximative zonotope
Z = zonotope(qZ);
ZmatFirst = get(Z,'Z');
dim = length(Q);

%order reduction after dependent generators
Zrest = zonotope([0*qZ.c,ZmatFirst(:,depGens+2:end)]);
load errorOrder
Zred = reduce(Zrest,'girard',errorOrder);

Zmat = ZmatFirst(:,1:depGens+1);
gens = length(Zmat(1,:)) - 1;

%for each dimension, compute generator elements
for i = 1:dim
    
    %pure quadratic evaluation
    quadMat = Zmat'*Q{i}*Zmat;
    
    %center, important: start from depGens+1(+1) 
    c(i,1) = quadMat(1,1);
    
    %generators with center
    ind = 1:gens;
    G(i,ind) = quadMat(1,ind+1) + quadMat(ind+1,1)';
    
    %generators from diagonal elements
    Gsquare(i,ind) = diag(quadMat(ind+1,ind+1));
    
    %generators from other elements
    counter = 0;
    
    for j = 1:gens
        kInd = (j+1):gens;
        Gquad(i, counter + kInd - j) = quadMat(j+1, kInd+1) + quadMat(kInd+1, j+1)';
        counter = counter + length(kInd);
    end
end

% quadratic multiplication of remaining generators
Grest = [];
Zquad_rest = quadraticMultiplication(Zred,Q);
G_add = generators(Zquad_rest);
Grest(:,end+1:end+length(G_add(1,:))) = G_add;

% mixed multiplications
Zrem = zonotope(Zmat);
Zquad_mixed1 = mixedMultiplication(Zred,Zrem,Q);
G_add = generators(Zquad_mixed1);
Grest(:,end+1:end+length(G_add(1,:))) = G_add;

Zquad_mixed2 = mixedMultiplication(Zrem,Zred,Q);
G_add = generators(Zquad_mixed2);
Grest(:,end+1:end+length(G_add(1,:))) = G_add;

%generate quadZonotope
qZquad = quadZonotope(c,G,Gsquare,Gquad,Grest);



%------------- END OF CODE --------------