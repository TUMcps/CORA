function M = supremum(intMat)
% supremum - returns the supremum of the interval matrix; mainly
%    implemented so that supremum can be used for both interval and
%    intervalMatrix objects without the need for case differentiation
%
% Syntax:
%    M = supremum(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    M - matrix representing the supremum
%
% Example:
%    intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
%    M = supremum(intMat);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix/infimum

% Authors:       Mark Wetzlinger
% Written:       16-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

M = supremum(intMat.int);

% ------------------------------ END OF CODE ------------------------------
