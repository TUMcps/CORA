function M = abs(intMat)
% abs - returns the absolute value bound of an interval matrix
%
% Syntax:
%    M = abs(intMat) 
%
% Inputs:
%    intMat - interval matrix
%
% Outputs:
%    M - absolute value bound
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Matthias Althoff
% Written:       21-July-2010
% Last update:   26-August-2011
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

M = supremum(abs(intMat.int));

% ------------------------------ END OF CODE ------------------------------
