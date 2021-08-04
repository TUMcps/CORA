function Z = box(Z)
% box - computes an enclosing axis-aligned box
%
% Syntax:  
%    Z = box(Z)
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
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:        Matthias Althoff
% Written:       09-March-2009
% Last update:   27-Aug-2019
% Last revision: ---

%------------- BEGIN CODE --------------

% determine new generator matrix
G = diag(sum(abs(generators(Z)),2));

% instantiate axis-aligned zonotope
Z.Z = [center(Z),G];

%------------- END OF CODE --------------