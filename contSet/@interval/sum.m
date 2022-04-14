function res = sum(obj,varargin)
% sum - Overloaded 'sum()' operator for intervals
%
% Syntax:  
%    res = sum(obj)
%    res = sum(obj,n)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      05-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
n = 1;
if nargin <= 1
    if size(obj,1) == 1
        n = 2;
    end
else
    n = varargin{1};
end

% init
res = interval();

% infimum
res.inf = sum(obj.inf,n);

% supremum
res.sup = sum(obj.sup,n);

%------------- END OF CODE --------------