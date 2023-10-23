function I = atan(I)
% atan - Overloaded 'atan()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [atan(x_), atan(x--)].
%
% Syntax:
%    I = atan(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I = interval([-1;1]);
%    res = atan(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       05-February-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                21-May-2022 (MW, remove new instantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I.inf = atan(I.inf);
I.sup = atan(I.sup);

% ------------------------------ END OF CODE ------------------------------
