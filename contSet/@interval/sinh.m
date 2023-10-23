function I = sinh(I)
% sinh - Overloaded 'sinh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [sinh(x_), sinh(x--)].
%
% Syntax:
%    I = sinh(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([-2;-3],[5;6]);
%    sinh(I)
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

I.inf = sinh(I.inf);
I.sup = sinh(I.sup);

% ------------------------------ END OF CODE ------------------------------
