function [Sigma]=sigma(pZ)
% sigma - returns Sigma matrix of a probabilistic zonotope
%
% Syntax:  
%    [Sigma]=sigma(pZ)
%
% Inputs:
%    pZ - probabilistic zonotope object
%
% Outputs:
%    Sigma - sigma matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 28-August-2007
% Last update: 26-February-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%reduce probabilistic zonotope first
[pZ]=probReduce(pZ);

%get new sigma matrix
G=pZ.g;
Sigma=G*G';

%------------- END OF CODE --------------