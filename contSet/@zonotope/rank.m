function n = rank(Z)
% rank - computes the rank of the generator matrix of a zonotope
%
% Syntax:  
%    n = rank(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    n - rank of the zonotope Z
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    n = rank(Z) % n=2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       06-May-2009
% Last update:   15-Sep-2019 (rename dim -> rank)
% Last revision: ---

%------------- BEGIN CODE --------------

n = rank(generators(Z));

%------------- END OF CODE --------------