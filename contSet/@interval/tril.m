function res = tril(I,varargin)
% tril - gets lower triangular part of I
%
% Syntax:
%    res = tril(I)
%    res = tril(I,K)
%
% Inputs:
%    I - interval object
%    K - (see built-in tril for matrices)
%
% Outputs:
%    res - lower triangular interval object
%
% Example: 
%    I = interval(-ones(2), ones(2));
%    tril(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       12-October-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(varargin)
    K = 0;
elseif length(varargin) == 1
    K = varargin{1};
    inputArgsCheck({{K,'att',{'double'},{'scalar','>=',0}}});
else
    throw(CORAerror('CORA:tooManyInputArgs',1));
end
res = I;
res.inf = tril(res.inf,K);
res.sup = tril(res.sup,K);

% ------------------------------ END OF CODE ------------------------------
