function val = infNorm(A)
% infNorm - returns the maximum of the infinity norm of an interval matrix
%
% Syntax:  
%    val = infNorm(A)
%
% Inputs:
%    A - interval matrix
%
% Outputs:
%    val - infinity norm of the interval matrix
%
% Example: 
%    A=interval(rand(3),rand(3)+1);
%    n=infNorm(A);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      12-February-2007 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

val=norm(abs(A),inf);

%------------- END OF CODE --------------