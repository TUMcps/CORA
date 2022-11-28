function res = ne(I1,I2)
% ne - overloads '~='-operator for intervals
%
% Syntax:  
%    res = ne(I1,I2)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%
% Outputs:
%    res - boolean 
%
% Example:
%    I1 = interval(-1,0);
%    I2 = interval(-2,-1);
%    I1 ~= I2
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Dmitry Grebenyuk
% Written:      06-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = any(infimum(I1) ~= infimum(I2)) ||...
    any(supremum(I1) ~= supremum(I2));

%------------- END OF CODE --------------