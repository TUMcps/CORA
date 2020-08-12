function d = dim(E)
% dim - returns the dimension of E
%
% Syntax:  
%    d = dim(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    d - dimension of E
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    d = dim(E) 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      15-Sep-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
d = E.dim;
%------------- END OF CODE --------------