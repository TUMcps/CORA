function n = dim(S)
% dim - returns the dimension of the ambient space of a continuous set;
%    currently only used for the empty set constructor contSet()
%
% Syntax:  
%    n = dim(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    S = contSet();
%    n = dim(S)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      17-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

n = S.dimension;

%------------- END OF CODE --------------