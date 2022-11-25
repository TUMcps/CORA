function [n]=infNorm(A)
% infNorm - returns the maximum of the infinity norm of an interval matrix
%
% Syntax:  
%    [n]=infNorm(A)
%
% Inputs:
%    A - interval matrix
%
% Outputs:
%    n - infinity norm of the interval matrix
%
% Example: 
%    A=interval(rand(3),rand(3)+1);
%    n=infNorm(A);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Matthias Althoff
% Written: 12-February-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

n=norm(abs(A),inf);

%------------- END OF CODE --------------