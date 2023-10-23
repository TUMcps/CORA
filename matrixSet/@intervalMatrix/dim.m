function n = dim(intMat,rc)
% dim - returns the dimension of the interval matrix
%
% Syntax:
%    n = dim(intMat)
%    n = dim(intMat,rc)
%
% Inputs:
%    intMat - intervalMatrix object
%    rc - 1 for row dimension, 2 for column dimension
%
% Outputs:
%    n - array with row and column dimension
%
% Example: 
%    intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
%    n = dim(intMat);
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

% choose wlog dimension of infimum
if nargin == 1
    n = size(intMat.int.inf);

else
    % read dimension of rows or columns
    n = size(intMat.int.inf,rc);

end

% ------------------------------ END OF CODE ------------------------------
