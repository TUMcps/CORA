function probZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a probabilistic zonotope according to [1,(4)]
%
% Syntax:
%    probZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical or interval matrix
%    factor2 - probabilistic zonotope
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

% Authors:       Matthias Althoff
% Written:       29-August-2007
% Last update:   27-September-2007
%                16-June-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Find a probabilistic zonotope object
[probZ,matrix] = findClassArg(factor1,factor2,'probZonotope');

try

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

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if representsa_(probZ,'emptySet',eps)
        return
    end

    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
