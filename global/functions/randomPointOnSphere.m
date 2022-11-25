function [x]=randomPointOnSphere(dim)
% randomVector - generates a random vector
%
% Syntax:  
%    [x]=randomVector(dim)
%
% Inputs:
%    dim - dimension
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

% Author: Matthias Althoff
% Written: 30-September-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%generate normally distributed random varibles
x = randn(dim,1);

%normalize result
x=x/norm(x);

%------------- END OF CODE --------------