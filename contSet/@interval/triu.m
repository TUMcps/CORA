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
res.inf = triu(res.inf,K);
res.sup = triu(res.sup,K);

% ------------------------------ END OF CODE ------------------------------
