function Z = abs(Z)
% abs - returns a zonotope with absolute values of the center and the
%    generators, i.e., Z = (|c|,|g_1|,...,|g_n|)
%
% Syntax:  
%    Z = abs(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z = zonotope([1 -1 0; 0 0 -1]);
%    Z = abs(Z); % Z=[1 1 0; 0 0 1]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       30-September-2006 
% Last update:   22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

Z.Z = abs(Z.Z);

%------------- END OF CODE --------------