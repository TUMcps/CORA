function res = isequal(h1,h2,tol)
% isequal - checks if halfspaces are equal
%
% Syntax:  
%    res = isequal(h1,h2,tol)
%
% Inputs:
%    h1 - halfspace object
%    h2 - halfspace object
%    tol - tolerance (optional)
%
% Outputs:
%    res - boolean whether h1 and h2 are equal
%
% Example: 
%    h1 = halfspace([2;4], 4);
%    h2 = halfspace([-3;5], 3);
%    isequal(h1,h2);
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

res = all(abs(h1.c - h2.c) < tol) && ... % normal vectors
    abs(h1.d - h2.d) < tol; % distances to origin
    

%------------- END OF CODE --------------