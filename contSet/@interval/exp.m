function I = exp(I)
% exp - Overloaded 'exp()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [exp(x_), exp(x--)].
%
% Syntax:
%    I = exp(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I = interval(-2,4);
%    exp(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   18-January-2024 (MW, avoid constructor call)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% exponential function is monotonic -> apply exp to infima/suprema
I.inf = exp(I.inf);
I.sup = exp(I.sup);

% ------------------------------ END OF CODE ------------------------------
