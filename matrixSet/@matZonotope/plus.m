function matZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of two matrix
%    zonotopes or a matrix zonotope with a matrix
%
% Syntax:
%    matZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - matZonotope object or numerical matrix
%    summand2 - matZonotope object or numerical matrix
%
% Outputs:
%    matZ - matrix zonotope after Minkowski addition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   05-August-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Find a matrix zonotope object
[matZ,summand] = findClassArg(summand1,summand2,'matZonotope');

%Is summand a zonotope?
if isa(summand,'matZonotope')
    %Calculate minkowski sum
    matZ.C = matZ.C + summand.center;
    matZ.G = cat(3,matZ.G,summand.G);
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    matZ.C = matZ.C + summand;
end

% ------------------------------ END OF CODE ------------------------------
