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
%    Z = zonotope([randn(2,1),randn(2,10)]);
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
