function intMat = intervalMatrix(matZ)
% intervalMatrix - computes an enclosing interval matrix of a matrix
%    zonotope
%
% Syntax:
%    intMat = intervalMatrix(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    intMat - intervalMatrix object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   06-October-2010
%                26-August-2011
%                03-April-2023 (MW, remove setting property)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% center matrix
C = center(matZ);

% delta matrix
D = abs(matZ.generator{1});
for i=2:matZ.gens
    D = D + abs(matZ.generator{i});
end

%instantiate interval matrix
intMat = intervalMatrix(C, D);


% ------------------------------ END OF CODE ------------------------------
