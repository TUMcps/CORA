function I = log(I)
% log - Overloaded log (natural logorithm) function for intervals 
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x-- < 0),
% [NaN, log(x--)] if (x_ < 0) and (x-- >= 0)
% [log(x_), log(x--)] if (x_ >= 0).
%
% Syntax:
%    I = log(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([3;9]);
%    res = log(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       07-February-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                05-May-2020 (MW, addition of error message)
%                21-May-2022 (MW, remove new instantiation)
%                18-January-2024 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute logarithm
I.inf = log(I.inf);
I.sup = log(I.sup);

% set entries with non-zero imaginary parts to NaN
I.inf(imag(I.inf) ~= 0) = NaN;
I.sup(imag(I.sup) ~= 0) = NaN;

% return error if NaN occures
if any(any(isnan(I.inf))) || any(any(isnan(I.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','>= 0'));
end

% ------------------------------ END OF CODE ------------------------------
