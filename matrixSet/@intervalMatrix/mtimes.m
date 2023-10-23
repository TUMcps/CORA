function intMat = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of matrix or an 
%    interval matrix with an interval matrix
%
% Syntax:
%    intMat = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or interval matrix
%    factor2 - numerical matrix or interval matrix
%
% Outputs:
%    intMat - interval matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   05-August-2010
%                03-April-2022 (MW, remove setting check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%factor1 is a numeric matrix
if isnumeric(factor1)
    %initialize factor
    matrix=factor1;
    %initialize matrix zonotope
    intMat=factor2;
    %compute new intervals
    intMat.int=matrix*intMat.int;
    
    
%factor2 is a numeric matrix
elseif isnumeric(factor2)
    %initialize factor
    matrix=factor2;
    %initialize matrix zonotope
    intMat=factor1;
    %compute new intervals
    intMat.int=intMat.int*matrix;
    
% multiplication with empty set
elseif isa(factor2,'emptySet')
    intMat = emptySet(dim(factor1,1));

% multiplication with full space
elseif isa(factor2,'fullspace')
    intMat = fullspace(dim(factor1,1));

%both factors are interval matrices
else
    %initialize result
    intMat=factor1;
    %compute interval matrix
    intMat.int=factor1.int*factor2.int;
end

% ------------------------------ END OF CODE ------------------------------
