function res = isempty(matP)
% isempty - checks if an interval matrix is empty
%
% Syntax:
%    res = isempty(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    res - true/false
%
% Example:
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2]; G{2} = [2 0; 1 -1];
%    matZ = matZonotope(C,G);
%    res = isempty(matZ);
%
%    matZ = matZonotope();
%    res = isempty(matZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check dimension of interval
res = any(dim(matP) == 0);

% ------------------------------ END OF CODE ------------------------------
