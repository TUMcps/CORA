function [res,S_conv] = representsa(matZ,varargin)
% representsa - checks if a set can also be represented by a different set,
%    e.g., a special case
%
% Syntax:
%    res = representsa(matZ,type)
%    res = representsa(matZ,type,tol)
%    [res,S_conv] = representsa(matZ,type)
%    [res,S_conv] = representsa(matZ,type,tol)
%
% Inputs:
%    matZ - matZonotope object
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

% reuse functionality from zonotope class
Z = zonotope(matZ);

if nargout == 1
    res = representsa(Z,varargin{:});
else
    [res,S_conv] = representsa(Z,varargin{:});
end

end

% ------------------------------ END OF CODE ------------------------------
