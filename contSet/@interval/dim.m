function d = dim(obj)
% dim - returns the dimension of the space the interval is a subset of
%
% Syntax:  
%    d = dim(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    d - dimension of obj 
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      15-Sep-2019
% Last update:  09-June-2020 (handling of interval matrices)
%               12-March-2021 (MW, empty intervals)
% Last revision:---

%------------- BEGIN CODE --------------

infi = infimum(obj); % equivalently: supremum(obj)
if isempty(infi)
    d = 0; return;
end

[rows, cols] = size(infi);
if rows == 1
    d = cols;
elseif cols == 1
    d = rows;
else
    d = [rows, cols];
end


%------------- END OF CODE --------------