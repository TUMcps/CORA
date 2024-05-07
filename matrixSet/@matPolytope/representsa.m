function [res,S_conv] = representsa(matP,varargin)
% representsa - checks if a set can also be represented by a different set,
%    e.g., a special case
%
% Syntax:
%    res = representsa(matP,type)
%    res = representsa(matP,type,tol)
%    [res,S_conv] = representsa(matP,type)
%    [res,S_conv] = representsa(matP,type,tol)
%
% Inputs:
%    matP - matPolytope object
%    type - char array
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%    S_conv - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       02-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reuse functionality from polytope class
P = polytope(matP);

if nargout == 1
    res = representsa(P,varargin{:});
else
    [res,S_conv] = representsa(P,varargin{:});
end

end

% ------------------------------ END OF CODE ------------------------------
