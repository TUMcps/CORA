function [res,S_conv] = representsa(intMat,varargin)
% representsa - checks if a set can also be represented by a different set,
%    e.g., a special case
%
% Syntax:
%    res = representsa(intMat,type)
%    res = representsa(intMat,type,tol)
%    [res,S_conv] = representsa(intMat,type)
%    [res,S_conv] = representsa(intMat,type,tol)
%
% Inputs:
%    S - contSet object
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
% Written:       26-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reuse functionality from interval class
I = interval(intMat);

% reshape to vector (required for respresentsa)
I = reshape(I,[],1);

if nargout == 1
    res = representsa(I,varargin{:});
else
    [res,S_conv] = representsa(I,varargin{:});
end

end

% ------------------------------ END OF CODE ------------------------------
