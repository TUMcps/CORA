function res = isempty(intMat)
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
%    intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
%    res = isempty(intMat);
%
%    intMat = intervalMatrix();
%    res = isempty(intMat);
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
res = any(dim(intMat) == 0);

% ------------------------------ END OF CODE ------------------------------
