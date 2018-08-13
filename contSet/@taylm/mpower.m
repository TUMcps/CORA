function res = mpower(base,exponent)
% mpower - Overloaded '^' operator for taylm (power)
%
% Syntax:  
%    res = mpower(base,exponent)
%
% Inputs:
%    base - taylm object
%    exponent - taylm object
%
% Outputs:
%    res - taylm
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      06-May-2016             
% Last update:  03-October-2017 (DG) A^0 = I is added               
% Last revision:---

%------------- BEGIN CODE --------------
if isscalar(base)
    res = power(base,exponent);
    
elseif exponent > 0 && all(size(base) == size(base'))
    res = base;
    for i = 2:exponent
        res = res * base;
    end
elseif exponent == 0 && all(size(base) == size(base'))
    res = eye(size(base,1));
else
    error('Matrix is not square')
end
end

