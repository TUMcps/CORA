function res = isequal(Int1,Int2,tol)
% isequal - checks if two intervals are equal
%
% Syntax:  
%    res = isequal(Int1,Int2,tol)
%
% Inputs:
%    Int1 - interval
%    Int2 - interval
%    tol - tolerance (optional)
%
% Outputs:
%    res - boolean whether intervals are equal
%
% Example: 
%    Int1 = interval([1; -1; 0], [4; 2; 1]);
%    Int2 = interval([1; 0; 0], [3.5; 2; 1]);
%    res = isequal(Int1,Int2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 2
    tol = eps;
end

res = all(abs(infimum(Int1) - infimum(Int2)) < tol) && ... % infima
    all(abs(supremum(Int1) - supremum(Int2)) < tol); % suprema


%------------- END OF CODE --------------