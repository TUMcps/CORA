function c = center(intMat)
% center - Returns the center of an intervalMatrix
%
% Syntax:
%    c = center(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    c - center of the interval matrix
%
% Example:
%    M = intervalMatrix(eye(2),2*eye(2));
%    c = center(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-May-2021 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = center(intMat.int);

% ------------------------------ END OF CODE ------------------------------
