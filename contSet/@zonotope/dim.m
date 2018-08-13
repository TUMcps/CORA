function [n] = dim(Z)
% dim - Returns the dimension of a zonotope in the sense that the rank of
% the generator matrix is computed
%
% Syntax:  
%    [n] = dim(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    n - dimension of the zonotope Z
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    n=dim(Z)-->n=2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-May-2009
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

n=rank(Z.Z(:,2:end));

%------------- END OF CODE --------------