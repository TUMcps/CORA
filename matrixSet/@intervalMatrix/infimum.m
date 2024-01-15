function M = infimum(intMat)
% infimum - returns the infimum of the interval matrix; mainly implemented
%    so that infimum can be used for both interval and intervalMatrix
%    objects without the need for case differentiation
%
% Syntax:
%    M = infimum(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    M - matrix representing the infimum
%
% Example:
%    intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
%    M = infimum(intMat);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

M = infimum(intMat.int);

% ------------------------------ END OF CODE ------------------------------
