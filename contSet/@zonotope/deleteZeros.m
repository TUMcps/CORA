function [Zred]=deleteZeros(Z)
% deleteZeros - removes zero generators
%
% Syntax:  
%    [Zred]=deleteZeros(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zred - reduced zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Matthias Althoff
% Written: 15-January-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get Z-matrix from zonotope Z
Zmatrix=get(Z,'Z');

%extract generator matrix
c=Zmatrix(:,1);
G=Zmatrix(:,2:end);

%Delete zero-generators
G=nonzeroFilter(G);

Zred=zonotope([c,G]);

%------------- END OF CODE --------------
