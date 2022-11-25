function res = mrdivide(numerator,denominator)
% mrdivide - Overload '/' operator
%
% Syntax:  
%    res = mrdivide(numerator,denominator)
%
% Inputs:
%    numerator - numerator (class zoo)
%    denominator - denominator (class zoo)
%
% Outputs:
%    res - resulting zoo object
%
% Other m-files required: inverse
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, interval

% Author:       Dmitry Grebenyuk
% Written:      06-November-2017
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE -------------

    res = numerator ./ denominator;
end

%------------ END OF CODE ------------