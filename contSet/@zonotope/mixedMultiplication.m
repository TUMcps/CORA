function [Zquad] = mixedMultiplication(Z1,Z2,Q)
% mixedMultiplication - computes \{Q_{ijk}*x_j*y_k|x \in Z1, y \in Z2\}
%
% Syntax:  
%    [Zquad] = mixedMultiplication(Z1,Z2,Q)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
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
% Written:      18-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%get matrix of zonotope
Zmat1 = get(Z1,'Z');
Zmat2 = get(Z2,'Z');
dim = length(Q);

%for each dimension, compute generator elements
for i = 1:dim
    
    %pure quadratic evaluation
    quadMat = Zmat1'*Q{i}*Zmat2;
    
    %vectorize
    quadVec = reshape(quadMat,1,[]);
    
    %center
    c(i,1) = quadVec(1);
    
    %generators
    G(i, :) = quadVec(2:end);
end

%generate new zonotope
Zquad = zonotope([c, G]);


%------------- END OF CODE --------------