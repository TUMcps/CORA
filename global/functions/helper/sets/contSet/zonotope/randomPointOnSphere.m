function x = randomPointOnSphere(n)
% randomPointOnSphere - generates a random vector
%
% Syntax:
%    x = randomPointOnSphere(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    x - random vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-September-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate normally distributed random variable
x = randn(n,1);

% normalize result
x = x/norm(x);

% ------------------------------ END OF CODE ------------------------------
