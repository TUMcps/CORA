function res = mrdivide(numerator,denominator)
% mrdivide - Overloaded matrix division '/' operator for intervals
%
% Syntax:
%    res = mrdivide(numerator, denominator)
%
% Inputs:
%    numerator - interval object
%    denominator - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    numerator = interval([-2;1],[3;2]);
%    denominator = 2;
%    numerator / denominator
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   01-July-2015
%                10-September-2015
%                13-March-2016 (speed improvement)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isscalar(denominator)
    res = numerator ./ denominator;
else
    throw(CORAerror('CORA:noops',numerator,denominator));
end

% ------------------------------ END OF CODE ------------------------------
