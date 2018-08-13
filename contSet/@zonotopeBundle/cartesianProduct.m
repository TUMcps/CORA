function Z1 = cartesianProduct(Z1,Z2)
% cartesianProduct - Returns the cartesian product of a zonotopeBundle with
%                    a zonotopes
%
% Syntax:  
%    Z = cartesianProduct(Z1,Z2)
%
% Inputs:
%    Z1 - zonotopeBundle object
%    Z2 - zonotope object
%
% Outputs:
%    Z1 - zonotopeBundle object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      13-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    for i = 1:Z1.parallelSets
       Z1.Z{i} = cartesianProduct(Z1.Z{i},Z2); 
    end

%------------- END OF CODE --------------