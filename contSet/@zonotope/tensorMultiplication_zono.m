function Zres = tensorMultiplication_zono(Z,M,options)
% tensorMultiplication_zono - computes \{M_{ijk...l}*x_j*x_k*...*x_l|x \in Z\}
%    when the center of Z is the origin and M is a matrix zonotope
%
% Syntax:
%    Zres = tensorMultiplication_zono(Z,M,options)
%
% Inputs:
%    Z - zonotope object
%    M - tensor
%    options - ?
%
% Outputs:
%    Zres - zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       10-October-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute zonotope of center
Zres = tensorMultiplication(Z,M.center,options);

%add results from generators
for i = 1:M.gens
    Zres = Zres + tensorMultiplication(Z,M.generator{i},options);
end

% ------------------------------ END OF CODE ------------------------------
