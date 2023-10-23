function res = prod(I,varargin)
% prod - product of array elements
%
% Syntax:
%    res = prod(I)
%    res = prod(I,n)
%
% Inputs:
%    I - interval object
%    n - dimension along which the product should be computed
%           1 - product of column's elements
%           2 - product of row's elements
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-3;-2;-4],[4;1;6]);
%    prod(I,1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Dmitry Grebenyuk
% Written:       24-October-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
n = 1;
if nargin <= 1
    if size(I,1) == 1
        n = 2;
    end
elseif ~isnumeric(varargin{1}) || ~isscalar(varargin{1}) || ~any(varargin{1} == [1,2])
    throw(CORAerror('CORA:wrongValue','second','either 1 or 2'));
else
    n = varargin{1};
end

% check input arguments
inputArgsCheck({{I,'att','interval'};
                {n,'att','numeric','nonnan'}});

if n == 1 % reduce to a row 
    S.type='()'; % to avoid Matlab's bug
    S.subs={1,':'};
    res = subsref(I,S);
    for i = 2:size(I, n)
        S.subs={i,':'};
        res = res .* subsref(I,S);
    end
elseif n == 2 % reduce to a column
    S.type='()';
    S.subs={':', 1};
    res = subsref(I,S);
    for i = 2:size(I, n)
        S.subs={':', i};
        res = res .* subsref(I,S);
    end
end

% ------------------------------ END OF CODE ------------------------------
