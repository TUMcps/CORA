function Z = cartesianProduct(Z1,Z2)
% cartesianProduct - Returns the cartesian product of two zonotopes
%
% Syntax:  
%    Z = cartesianProduct(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      18-May-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c = [Z1.Z(:,1); Z2.Z(:,1)];

%determine sizes of generator matrices
[rows_1, cols_1] = size(Z1.Z(:,2:end));
[rows_2, cols_2] = size(Z2.Z(:,2:end));

%copy generators
G(1:rows_1, 1:cols_1) = Z1.Z(:,2:end);
G((rows_1+1):(rows_1+rows_2), (cols_1+1):(cols_1+cols_2)) = Z2.Z(:,2:end);

%generate new zonotope
Z = zonotope([c,G]);

%------------- END OF CODE --------------