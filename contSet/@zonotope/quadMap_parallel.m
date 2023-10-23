function Z = quadMap_parallel(Z,Q)
% quadMap_parallel - computes \{Q_{ijk}*x_j*x_k|x \in Z\} using parfor
%
% Syntax:
%    Z = quadMap_parallel(Z,Q)
%
% Inputs:
%    Z - zonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z = zonotope([0 1 1;0 1 0]);
%    Q{1} = [0.5 0.5; 0 -0.5];
%    Q{2} = [-1 0; 1 1];
% 
%    res = quadMap_parallel(Z,Q);
% 
%    figure; hold on;
%    plot(Z,[1,2],'r');
% 
%    figure; hold on;
%    plot(res,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-December-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get center and generator matrix of zonotope
Zmat = [Z.c,Z.G];
dimQ = length(Q);
gens = size(Zmat,2) - 1;

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

% generate new zonotope
Z = zonotope(c, Gmat);

% ------------------------------ END OF CODE ------------------------------
