function c = center(obj)
% center - Returns the center of an intervalMatrix
%
% Syntax:  
%    c = center(obj)
%
% Inputs:
%    obj - intervalMatrix object
%
% Outputs:
%    c - center of the intervalMatrix obj
%
% Example:
%    M = intervalMatrix(eye(2),2*eye(2));
%    c = center(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-May-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

c = center(obj.int);

%------------- END OF CODE --------------