function [qZ] = quadZonotope(Z)
% quadZonotope - converts a zonotope to a quadZonotope format
%
% Syntax:  
%    [qZ] = quadZonotope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    qZ - quadZonotope object
%
% Example: 
%    ---
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      05-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%extract center and generators
c=Z.Z(:,1);
G=Z.Z(:,2:end);

%instantiate quadZonotope
qZ=quadZonotope(c,G,[],[],[]);

%------------- END OF CODE --------------