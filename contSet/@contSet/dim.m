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

% Authors:       Mark Wetzlinger
% Written:       17-June-2022
% Last update:   22-March-2023 (MW, adapt to new constructor syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = 0;

% ------------------------------ END OF CODE ------------------------------
