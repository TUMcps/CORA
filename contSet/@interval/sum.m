function res = sum(I,varargin)
% sum - Overloaded 'sum()' operator for intervals
%
% Syntax:
%    res = sum(I)
%    res = sum(I,n)
%
% Inputs:
%    I - interval object
%    n - dimension along which the sum should be computed
%           1 - sum of column's elements
%           2 - sum of row's elements
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-3;-2;-4],[4;1;6]);
%    sum(I,1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-August-2016
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

% init (overwritten below)
res = interval.empty(1);

% infimum
res.inf = sum(I.inf,n);

% supremum
res.sup = sum(I.sup,n);

% ------------------------------ END OF CODE ------------------------------
