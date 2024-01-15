function n = dim(I)
% dim - returns the dimension of the ambient space of an interval
%
% Syntax:
%    n = dim(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example:
%    I = interval([-1;-2],[3;4]);
%    n = dim(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       15-September-2019
% Last update:   09-June-2020 (handling of interval matrices)
%                12-March-2021 (MW, empty intervals)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

infi = infimum(I); % equivalently: supremum(I)
if isempty(infi)
    n = size(infi,1); return;
end

[rows, cols] = size(infi);
if rows == 1
    n = cols;
elseif cols == 1
    n = rows;
else
    n = [rows, cols];
end


% ------------------------------ END OF CODE ------------------------------
