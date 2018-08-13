function [Z] = mean(pZ)
% mean - Returns the uncertain mean of a probabilistic zonotope
%
% Syntax:  
%    [Z] = mean(pZ)
%
% Inputs:
%    pZ - probabilistic zonotope object
%
% Outputs:
%    Z - uncertain mean of the zonotope pZ
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 27-September-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

Z=zonotope(pZ.Z);

%------------- END OF CODE --------------