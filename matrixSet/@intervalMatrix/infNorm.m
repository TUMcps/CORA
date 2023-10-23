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
%    A = intervalMatrix([0 1; 1 -2],[1 1; 2 0]);
%    val = infNorm(A);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-February-2007 
% Last update:   13-April-2023 (MW, move from linProbSys/private -> intervalMatrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

val = norm(abs(A),Inf);

% ------------------------------ END OF CODE ------------------------------
