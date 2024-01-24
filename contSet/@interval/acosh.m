function I = acosh(I)
% acosh - Overloaded 'acosh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x-- < 1),
% [NaN, acosh(x--)] if (x_ < 1) and (x-- >= 1),
% [acosh(x_), acosh(x--)] if (x_ >= 1).
%
% Syntax:
%    I = acosh(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([2;4]);
%    res = acosh(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       12-February-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                05-May-2020 (MW, standardized error message)
%                21-May-2022 (MW, remove new instantiation)
%                18-January-2024 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% to preserve the shape
lb = I.inf;
ub = I.sup;

% find indices
ind1 = lb >= 1;
I.inf(ind1) = acosh(lb(ind1));
I.sup(ind1) = acosh(ub(ind1));

ind2 = lb < 1 & ub >= 1;
I.inf(ind2) = NaN;
I.sup(ind2) = acosh(I.sup(ind2));

ind3 = ub < 1;
I.inf(ind3) = NaN;
I.sup(ind3) = NaN;

% return error if NaN occures
if any(any(isnan(I.inf))) || any(any(isnan(I.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','>= -1 && <= 1'));
end

% ------------------------------ END OF CODE ------------------------------
