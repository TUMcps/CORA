function dp = fpolyder(p)
% fpolyder - computes the derivative of a given polynomial
%    same as polyder, but much faster
%
% Syntax:
%    dp = fpolyder(p)
%
% Inputs:
%    p - coefficient of polynomial as row vector
%
% Outputs:
%    dp - coefficient of derivative
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyder

% Authors:       Tobias Ladner
% Written:       30-May-2023
% Last update:   02-October-2023 (rename to fpolyder)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

order = length(p)-1;
dp = (order:-1:1) .* p(1:end-1);

% ------------------------------ END OF CODE ------------------------------
