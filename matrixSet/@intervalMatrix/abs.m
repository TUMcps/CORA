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
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author:       Matthias Althoff
% Written:      21-July-2010
% Last update:  26-August-2011
% Last revision:---

%------------- BEGIN CODE --------------

Minf=infimum(intMat.int);
Msup=supremum(intMat.int);

M=max(abs(Minf), abs(Msup));

%------------- END OF CODE --------------