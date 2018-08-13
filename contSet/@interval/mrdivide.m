function res = mrdivide(numerator,denominator)
% mrdivide - Overloaded matrix division '/' operator for intervals
%
% Syntax:  
%    res = mrdivide(numerator, denominator)
%
% Inputs:
%    numerator, denominator - interval objects
%
% Outputs:
%    res - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  01-July-2015
%               10-September-2015
%               13-March-2016           Speed improvement
% Last revision:---

%------------- BEGIN CODE --------------

% to initiate res
res = interval();

if isscalar(denominator)
    res = numerator ./ denominator;
else
    error('The function of mrdivide works only for an interval / a number case.')
end

%------------- END OF CODE --------------