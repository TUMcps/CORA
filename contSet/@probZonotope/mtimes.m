function probZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a probabilistic zonotope according to [1,(4)]
%
% Syntax:  
%    probZ = mtimes(matrix,probZ)
%
% Inputs:
%    matrix - numerical or interval matrix
%    probZ - probabilistic zonotope
%
% Outputs:
%    probZ - probZonotope after multiplication of a matrix with a zonotope
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    matrix = rand(2,2);
%    probZ * matrix;
%
% References:
%    [1] M. Althoff et al. "Safety assessment for stochastic linear systems 
%        using enclosing hulls of probability density functions", ECC 2009
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      29-August-2007
% Last update:  27-September-2007
%               16-June-2016
% Last revision: ---

%------------- BEGIN CODE --------------

%Find a probabilistic zonotope object
%Is factor1 a probabilistic zonotope?
if isa(factor1,'probZonotope')
    %initialize resulting probabilistic zonotope
    probZ=factor1;
    %initialize other summand
    matrix=factor2;
%Is factor2 a probabilistic zonotope?    
elseif isa(factor2,'probZonotope')
    %initialize resulting probabilistic zonotope
    probZ=factor2;
    %initialize other summand
    matrix=factor1;  
end

%numeric matrix
if isnumeric(matrix)
    probZ.Z=matrix*probZ.Z;
    %pZ.g=matrix*pZ.g;
    probZ.cov=matrix*probZ.cov*matrix';
    
%interval matrix
elseif isa(matrix,'interval')
    %get center of interval matrix
    T=center(matrix);
    
    %get symmetric interval matrix
    M_min = infimum(matrix);
    M_max = supremum(matrix);
    S=0.5*(M_max-M_min);
    
    %probabilistic zonotope to zonotope
    mSigmaZ=zonotope(probZ);
    Z=mSigmaZ.Z;
    Zsum=sum(abs(Z),2);    
    
    %compute new zonotope
    probZ.Z=[T*probZ.Z,diag(S*Zsum)]; 
    probZ.g=T*probZ.g;
end

%------------- END OF CODE --------------
