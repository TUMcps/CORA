function O = projectHighDim_(O,N,proj)
% projectHighDim_ - lifts an empty set onto a higher-dimensional space
%
% Syntax:
%    O = projectHighDim_(O,N,dims)
%
% Inputs:
%    O - emptySet object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    O - lifted emptySet object
%
% Example: 
%    O = emptySet(4);
%    projectHighDim(O,6,[1,2,5,6])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% function expects higher dimension to be bounded at 0, however, the set
% will be empty as it is empty in the existing dimensions
O.dimension = N;

% ------------------------------ END OF CODE ------------------------------
