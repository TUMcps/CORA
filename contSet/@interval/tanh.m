function I = tanh(I)
% tanh - Overloaded 'tanh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [tanh(x_), tanh(x--)].
%
% Syntax:
%    I = tanh(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([-2;-3],[3;4]);
%    tanh(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       05-February-2016
% Last update:   22-February-2016 (DG, the matrix case is rewritten)
%                21-May-2022 (MW, remove new instantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I.inf = tanh(I.inf);
I.sup = tanh(I.sup);

% ------------------------------ END OF CODE ------------------------------
