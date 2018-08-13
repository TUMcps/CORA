function [Z] = box(Z)
% box - Computes an enclosing axis-aligned box
%
% Syntax:  
%    [Z] = box(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z=zonotope(rand(2,5));
%    B=box(Z);
%    plot(Z);
%    hold on
%    plot(B);
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author: Matthias Althoff
% Written: 09-March-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%extract generators
G=Z.Z(:,2:end);

%determine new generator matrix
G=diag(sum(abs(G),2));

%instantiate interval hull
Z.Z=[Z.Z(:,1),G];

%------------- END OF CODE --------------