function I = cosh(I)
% cosh - Overloaded 'cosh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [cosh(x--), cosh(x_)] if (x-- <= 0),
% [cosh(x_), cosh(x--)] if (x_ >= 0),
% [1, max( acosh(x_), acosh(x--)] if (x_ < 0) and (x-- > 0).
%
% Syntax:
%    I = cosh(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I = interval([-2;3],[3;4]);
%    res = cosh(I);
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
%                18-January-2024 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

lb = I.inf;
ub = I.sup;

% case 1: upper bound smaller than 0
ind1 = ub <= 0;
I.inf(ind1) = cosh(ub(ind1));
I.sup(ind1) = cosh(lb(ind1));

% case 1: lower bound larger than 0
ind2 = lb >= 0;
I.inf(ind2) = cosh(lb(ind2));
I.sup(ind2) = cosh(ub(ind2));

% remaining entries
ind3 = ~ind1 & ~ind2;
I.inf(ind3) = 1;
I.sup(ind3) = cosh( max( abs(lb(ind3)), abs(ub(ind3)) ) );

% ------------------------------ END OF CODE ------------------------------
