function matP = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or a 
%    matrix polytope with a matrix polytope
%
% Syntax:
%    matP = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or matrix polytope
%    factor2 - numerical matrix or matrix polytope
%
% Outputs:
%    matP - matrix polytope
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       21-June-2010 
% Last update:   02-May-2024 (TL, new structure of V)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%factor1 is a numeric matrix
if isnumeric(factor1)
    %initialize factor
    matrix=factor1;
    %initialize matrix polytope
    matP=factor2;
    %compute vertices
    matP.V = pagemtimes(matrix, matP.V);
    
%factor2 is a numeric matrix
elseif isnumeric(factor2)
    %initialize factor
    matrix=factor2;
    %initialize matrix polytope
    matP=factor1;
    %compute vertices
    matP.V = pagemtimes(matP.V,matrix);
    
%both factors are polytope matrices
else
    % get vertices of first matPolytope
    matP1=factor1;
    V1 = matP1.V;
    [n1,m1,h1] = size(V1);
    % get vertices of second matPolytope
    matP2=factor2;
    V2 = matP2.V;
    [n2,m2,h2] = size(V2);
    
    % reshape Z2
    V2 = reshape(V2,n2,m2,1,h2);

    % multiply each matrix from either set
    V = pagemtimes(V1,V2);

    % reshape back ---

    % fix dimensions for scalar multiplication
    if n1 == 1 && m1 == 1
        n1 = n2;
    end
    if n2 == 1 && m2 == 1
        m2 = m1;
    end
    
    %  reshape
    V = reshape(V,n1,m2,h1*h2);

    % init resulting matrix polytope
    matP = matPolytope(V);
end

% ------------------------------ END OF CODE ------------------------------
