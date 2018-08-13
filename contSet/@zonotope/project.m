function [Z] = project(Z,dim)
% project - Returns a zonotope which is projected onto the specified
% dimensions
%
% Syntax:  
%    [Z] = project(Z,dim)
%
% Inputs:
%    Z - zonotope object
%    dim - projected dimensions
%
% Outputs:
%    Z - zonotope, whereas Z=(|c|,|g_1|,...,|g_n|)
%
% Example: 
%    Z=zonotope([1 -1 0; 0 0 -1]);
%    Z=abs(Z)-->Z=[1 1 0; 0 0 1]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 15-September-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

Z.Z=Z.Z(dim,:);

%------------- END OF CODE --------------