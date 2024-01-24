function G = nonzeroFilter(G,varargin)
% nonzeroFilter - filters out generators of length 0
%
% Syntax:
%    G = nonzeroFilter(G)
%    G = nonzeroFilter(G,tol)
%
% Inputs:
%    G - matrix of generators
%    tol - (optional) tolerance
%
% Outputs:
%    G - reduced matrix of generators
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       12-September-2008
% Last update:   02-September-2019
%                24-January-2024 (MW, add tolerance)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% delete zero-generators
G = G(:,any(G,1));

if nargin > 1
    % read out like this for speed reasons
    tol = varargin{1};
    G = G(:,vecnorm(G,2,1) > tol);
end

% ------------------------------ END OF CODE ------------------------------
