function res = sum(obj,varargin)
% sum - sum of array elements
%
% Syntax:
%    res = sum(obj)
%    res = sum(obj,n)
%
% Inputs:
%    obj - taylm object
%    n - dimension:
%           1 - sum of column's elements
%           2 - sum of row's elements
%
% Outputs:
%    res - resulting taylm object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/sum

% Authors:       Niklas Kochdumper
% Written:       20-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
n = 1;
if nargin <= 1
    if size(obj,1) == 1
        n = 2;
    end
else
    n = varargin{1};
end

inputArgsCheck({{obj,'att','taylm'};
                {n,'att','numeric','nonnan'}});

if n == 1 % reduce to a row 
    S.type='()'; % to avoid Matlab's bug
    S.subs={1,':'};
    res = subsref(obj,S);
    for i = 2:size(obj, n)
        S.subs={i,':'};
        res = res + subsref(obj,S);
    end
elseif n == 2 % reduce to a column
    S.type='()';
    S.subs={':', 1};
    res = subsref(obj,S);
    for i = 2:size(obj, n)
        S.subs={':', i};
        res = res + subsref(obj,S);
    end
else
    error ('Wrong input')
end

% ------------------------------ END OF CODE ------------------------------
