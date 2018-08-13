function [pZ] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
% interval matrix with a zonotope
%
% Syntax:  
%    [Z] = mtimes(matrix,Z)
%
% Inputs:
%    matrix - numerical or interval matrix
%    Z - zonotope object 
%
% Outputs:
%    Z - Zonotpe after multiplication of a matrix with a zonotope
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    matrix=[0 1; 1 0];
%    plot(Z);
%    hold on
%    Z=matrix*Z;
%    plot(Z);
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
if strcmp('probZonotope',class(factor1))
    %initialize resulting probabilistic zonotope
    pZ=factor1;
    %initialize other summand
    matrix=factor2;
%Is factor2 a probabilistic zonotope?    
elseif strcmp('probZonotope',class(factor2))
    %initialize resulting probabilistic zonotope
    pZ=factor2;
    %initialize other summand
    matrix=factor1;  
end

%numeric matrix
if isnumeric(matrix)
    pZ.Z=matrix*pZ.Z;
    %pZ.g=matrix*pZ.g;
    pZ.cov=matrix*pZ.cov*matrix';
    
%interval matrix
elseif strcmp('interval',class(matrix))
    %get center of interval matrix
    T=mid(matrix);
    
    %get symmetric interval matrix
    M_min = infimum(matrix);
    M_max = supremum(matrix);
    S=0.5*(M_max-M_min);
    
    %probabilistic zonotope to zonotope
    mSigmaZ=zonotope(pZ);
    Z=get(mSigmaZ,'Z');
    Zsum=sum(abs(Z),2);    
    
    %compute new zonotope
    pZ.Z=[T*pZ.Z,diag(S*Zsum)]; 
    pZ.g=[T*pZ.g];
end

%------------- END OF CODE --------------
