function res = rad(intMat)
% rad - returns the radius of an intervalMatrix
%
% Syntax:
%    res = rad(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    res - numerical value (matrix)
%
% Example: 
%    M = intervalMatrix(eye(2), 2*eye(2));
%    b = rad(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       06-May-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = rad(intMat.int);

% ------------------------------ END OF CODE ------------------------------
