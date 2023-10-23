function n = dim(tay)
% dim - returns the dimension of the ambient space of a taylor model
%
% Syntax:
%    n = dim(tay)
%
% Inputs:
%    tay - taylm object
%
% Outputs:
%    n - dimension
%
% Example: 
%    tay = taylm(interval(2,3));
%    n = dim(tay);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       02-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(tay)
    n = 0; return;
end

[rows, cols] = size(tay);
if rows == 1
    n = cols;
elseif cols == 1
    n = rows;
else
    n = [rows, cols];
end

% ------------------------------ END OF CODE ------------------------------
