function res = diag(obj)
% diag - Create diagonal matrix or get diagonal elements of matrix
%
% Syntax:  
%    res = diag(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - diagonal matrix or diagonal elements of matrix
%
% Example: 
%    I = interval([1 2; -1 1], [2 3; 0 10]);
%    d = diag(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Matthias Althoff
% Written:      02-November-2017 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain result
res = interval(diag(obj.inf), diag(obj.sup));

%------------- END OF CODE --------------