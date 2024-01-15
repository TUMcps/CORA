function c = center(I)
% center - returns the center of an interval
%
% Syntax:
%    c = center(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    c - center of interval (vector)
%
% Example: 
%    I = interval([-1;1],[1;2]);
%    c = center(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       26-June-2015
% Last update:   02-September-2019 (rename mid -> center)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set check
if representsa_(I,'emptySet',eps)
    c = zeros(dim(I),0); return ;
end

% compute center
c = 0.5*(I.inf + I.sup);

% for unbounded dimensions, return NaN
c(isinf(c)) = NaN;

% ------------------------------ END OF CODE ------------------------------
