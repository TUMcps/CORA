function O = lift_(O,N,dims)
% lift_ - lifts an empty set onto a higher-dimensional space
%
% Syntax:
%    O = lift_(O,N,dims)
%
% Inputs:
%    O - emptySet object
%    N - dimension of the higher-dimensional space
%    dims - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    O - lifted emptySet object
%
% Example: 
%    O = emptySet(4);
%    lift(O,6,[1,2,5,6])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift

% Authors:       Mark Wetzlinger
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% function expects higher dimension to be unbounded, however, the set
% will be empty as it is empty in the existing dimensions
O.dimension = N;

% ------------------------------ END OF CODE ------------------------------
