function res = triu(I,varargin)
% triu - gets upper triangular part of I
%
% Syntax:
%    res = triu(I)
%    res = triu(I,K)
%
% Inputs:
%    I - interval object
%    K - (see built-in tril for matrices)
%
% Outputs:
%    res - upper triangular interval object
%
% Example: 
%    I = interval(-ones(2), ones(2));
%    triu(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       12-October-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,2);

if nargin == 1
    K = 0;
elseif nargin == 2
    K = varargin{1};
    inputArgsCheck({{K,'att',{'double'},{'scalar',@(K) K>=0}}});
end

res = I;
res.inf = triu(res.inf,K);
res.sup = triu(res.sup,K);

% ------------------------------ END OF CODE ------------------------------
