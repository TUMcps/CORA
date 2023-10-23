function I = asinh(I)
% asinh - Overloaded 'asinh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [asinh(x_), asinh(x--)].
%
% Syntax:
%    I = asinh(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I = interval([-1;1]);
%    res = asinh(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       12-February-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                21-May-2022 (MW, remove new instantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I.inf = asinh(I.inf);
I.sup = asinh(I.sup);

% ------------------------------ END OF CODE ------------------------------
