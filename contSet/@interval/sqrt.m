function I = sqrt(I)
% sqrt - Overloaded 'sqrt'-function for intervals, computes the square root
%    of the interval
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x-- < 0),
% [NaN, sqrt(x--)] if (x_ < 0) and (x-- >= 0),
% [sqrt(x_), sqrt(x--)] if (x_ >= 0).
%
% Syntax:
%    res = sqrt(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval(9,16);
%    sqrt(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       20-January-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                05-May-2020 (MW, standardized error message)
%                18-January-2024 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% to preserve the shape
lb = I.inf;
ub = I.sup;

% find indices
ind1 = lb >= 0;
I.inf(ind1) = sqrt(lb(ind1));
I.sup(ind1) = sqrt(ub(ind1));

ind2 = lb < 0 & ub >= 0;
I.inf(ind2) = NaN;
I.sup(ind2) = sqrt(ub(ind2));

ind3 = ub < 0;
I.inf(ind3) = NaN;
I.sup(ind3) = NaN;

% return error if NaN occures
if any(any(isnan(I.inf))) || any(any(isnan(I.sup)))
	throw(CORAerror('CORA:outOfDomain','validDomain','>= 0'));
end

% ------------------------------ END OF CODE ------------------------------
