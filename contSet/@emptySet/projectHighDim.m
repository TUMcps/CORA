function O = projectHighDim(O,N,dims)
% projectHighDim - projects an empty set onto a higher-dimensional space
%
% Syntax:  
%    O = project(O,N,dims)
%
% Inputs:
%    O - emptySet object
%    N - dimension of the higher-dimensional space
%    dims - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    O - projected emptySet object
%
% Example: 
%    O = emptySet(4);
%    projectHighDim(O,6,[1,2,5,6])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      06-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check input arguments
inputArgsCheck({{O,'att','emptySet'};
                {N,'att','numeric',{'scalar','nonnegative','integer','>=',O.dimension}};
                {dims,'att','numeric',{'vector','nonnegative','integer'}}});

O.dimension = N;

%------------- END OF CODE --------------