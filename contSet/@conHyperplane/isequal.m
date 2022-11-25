function res = isequal(hyp1,hyp2)
% isequal - checks if conHyperplane objects are equal
%
% Syntax:  
%    res = isequal(h)
%
% Inputs:
%    hyp1 - conHyperplane object
%    hyp2 - conHyperplane object
%    tol  - tolerance (optional)
%
% Outputs:
%    res - boolean whether hyp is empty or not
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
% Written:      16-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 2
    tol = eps;
end

res = isequal(hyp1.h,hyp2.h,tol) && ... % halfspaces
    all(abs(hyp1.C - hyp2.C) < tol) && ... % C matrices
    abs(hyp1.d - hyp2.d) < tol; % distances

%------------- END OF CODE --------------
