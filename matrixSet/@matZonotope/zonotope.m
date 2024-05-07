function Z = zonotope(matZ)
% zonotope - Converts a matrix zonotope to a zonotope 
%
% Syntax:
%    Z = zonotope(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       18-June-2010 
% Last update:   25-April-2024 (TL, much faster implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get center
c = reshape(matZ.C,[],1);

% get generators
[m,n,h] = size(matZ.G);
G = reshape(matZ.G,n*m,h);
    
%instantiate zonotope
Z=zonotope([c,G]);

% ------------------------------ END OF CODE ------------------------------
