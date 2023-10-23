function str = string(I)
% string - conversion of an interval to a string
%
% Syntax:
%    str = string(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    str - string
%
% Example:
%    I = interval([-1;-2],[2;1]);
%    str = string(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out rows and columns
[rows,cols] = size(I);

str = strings([rows,cols]);
for i=1:rows
    for j=1:cols
        str(i,j) = "[" + I.inf(i,j) + "," + I.sup(i,j) + "]";
    end
end

% ------------------------------ END OF CODE ------------------------------
