function matZ = matZonotope(Z)
% matZonotope - converts the given zonotope to a matrix zonotope
%
% Syntax:
%    matZ = matZonotope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    matZ - matZonotope object
%
% Example: 
%    Z = zonotope([ 0.578 ; 0.433 ], [ 0.768 -0.642 0.248 0.606 0.962 -0.536 0.215 -0.185 0.096 -0.583 ; -0.214 0.267 -0.344 0.999 -0.746 -0.953 -0.778 0.768 -0.262 -0.118 ]);
%    matZ = matZonotope(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: matZonotope

% Authors:       Tobias Ladner
% Written:       25-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract zonotope properties
c = center(Z);
G = generators(Z);

% extend generators
[n,h] = size(G);
G = reshape(G,n,1,h);

% construct matrix zonotope
matZ = matZonotope(c,G);

end

% ------------------------------ END OF CODE ------------------------------
