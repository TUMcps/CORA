function [c] = center(Z)
% center - Returns the center of a zonotope
%
% Syntax:  
%    [c] = center(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    c - center of the zonotope Z
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    c=center(Z)-->c=[1;0]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 30-September-2006 
% Last update: 22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

c=Z.Z(:,1);

%------------- END OF CODE --------------