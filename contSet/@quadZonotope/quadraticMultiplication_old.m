function [qZquad] = quadraticMultiplication_old(qZ,Q)
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
% Last update:  ---
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
ZmatRest = get(Zred,'Z');

Zmat = [ZmatFirst(:,1:depGens), ZmatRest];
gens = length(Zmat(1,:)) - 1;

%for each dimension, compute generator elements
for i = 1:dim
    
    %pure quadratic evaluation
    quadMat = Zmat'*Q{i}*Zmat;
    
    %center, important: start from depGens+1(+1) 
    c(i,1) = quadMat(1,1) + 0.5*sum(diag(quadMat(depGens+2 : end,depGens+2 : end)));
    
    %generators with center
    indDep = 1:depGens;
    indInd = depGens+1 : gens;
    G(i,indDep) = quadMat(1,indDep+1) + quadMat(indDep+1,1)'; 
    Grest(i,indInd - depGens) = quadMat(1,indInd+1) + quadMat(indInd+1,1)'; 
    
    %generators from diagonal elements
    Gsquare(i,indDep) = diag(quadMat(indDep+1,indDep+1));
    Grest(i,gens + indInd - 2*depGens) = 0.5*diag(quadMat(indInd+1,indInd+1));
    
    %generators from other elements
    counter = 0;
    counterRest = 2*(gens-depGens);
    
    %1:depGens
    for j = 1:depGens
        kInd = (j+1) : depGens;
        kIndRest = (depGens+1) : gens;
        Gquad(i, counter + kInd - j) = quadMat(j+1, kInd+1) + quadMat(kInd+1, j+1)';
        Grest(i, counterRest + kIndRest - depGens) = quadMat(j+1, kIndRest+1) + quadMat(kIndRest+1, j+1)';
        counter = counter + length(kInd);
        counterRest = counterRest + length(kIndRest);
    end 
    
    %depGens+1:gens
    for j = (depGens+1):gens 
        kInd = (j+1):gens;
        Grest(i, counterRest + kInd - j) = quadMat(j+1, kInd+1) + quadMat(kInd+1, j+1)';
        counterRest = counterRest + length(kInd);
    end 
    
end

%generate quadZonotope
qZquad = quadZonotope(c,G,Gsquare,Gquad,Grest);



%------------- END OF CODE --------------