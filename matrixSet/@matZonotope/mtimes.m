function matZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or a 
%    matrix zonotope with a matrix zonotope
%
% Syntax:
%    matZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or matZonotope object
%    factor2 - numerical matrix or matZonotope object
%
% Outputs:
%    matZ - matrix zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       18-June-2010 
% Last update:   05-August-2010
%                25-April-2024 (TL, pagemtimes, much faster implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%factor1 is a numeric matrix
if isnumeric(factor1)
    %initialize factor
    matrix=factor1;
    %initialize matrix zonotope
    matZ=factor2;
    %compute center
    matZ.C=matrix*matZ.C;
    matZ.G = pagemtimes(matrix,matZ.G);
    
%factor2 is a numeric matrix
elseif isnumeric(factor2)
    %initialize factor
    matrix=factor2;
    %initialize matrix zonotope
    matZ=factor1;
    %compute center
    matZ.C=matZ.C*matrix;
    matZ.G = pagemtimes(matZ.G,matrix);
    
%both factors are zonotope matrices
else
    % initialize matrix zonotope 1
    matZ1=factor1;
    % initialize matrix zonotope 2
    matZ2=factor2;

    % concat center with generators
    Z1 = cat(3,matZ1.C,matZ1.G);
    [n1,m1,h1] = size(Z1);
    Z2 = cat(3,matZ2.C,matZ2.G);
    [n2,m2,h2] = size(Z2);

    % reshape Z2
    Z2 = reshape(Z2,n2,m2,1,h2);

    % multiply each matrix from either set
    Z = pagemtimes(Z1,Z2);

    % reshape back ---

    % fix dimensions for scalar multiplication
    if n1 == 1 && m1 == 1
        n1 = n2;
    end
    if n2 == 1 && m2 == 1
        m2 = m1;
    end
    
    %  reshape
    Z = reshape(Z,n1,m2,h1*h2);

    % construct final matrix zonotope
    matZ = matZonotope();
    matZ.C = Z(:,:,1);
    matZ.G = Z(:,:,2:end);
end

% ------------------------------ END OF CODE ------------------------------
 