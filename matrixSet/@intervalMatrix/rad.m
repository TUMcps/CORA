function res = rad(obj)
% rad - returns the radius of an intervalMatrix
%
% Syntax:  
%    res = rad(obj)
%
% Inputs:
%    obj - intervalMatrix object
%
% Outputs:
%    res - numerical value (matrix)
%
% Example: 
%    M = intervalMatrix(eye(2), 2*eye(2));
%    b = rad(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Victor Gassmann
% Written:      23-July-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0.5*(obj.Sup - obj.Inf);

%------------- END OF CODE --------------