function res = isequal(C1,C2,tol)
% isequal - checks if two capsules are equal
%
% Syntax:  
%    res = isequal(C1,C2)
%
% Inputs:
%    C1 - capsule
%    C2 - capsule
%    tol - tolerance (optional)
%
% Outputs:
%    res - boolean whether capsules are equal
%
% Example: 
%    C1 = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    C2 = capsule([1; 0; 0], [0.5; -1; 1], 0.5);
%    res = isequal(C1,C2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  12-March-2021 (MW, add dimension mismatch)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 2
    tol = eps;
end

if dim(C1) ~= dim(C2)
    [id,msg] = errDimMismatch();
    error(id,msg);
end

res = all(abs(center(C1) - center(C2)) < tol) && ... % center
    all(abs(C1.g - C2.g) < tol) && ... % generator
    abs(C1.r - C2.r) < tol; % radius


%------------- END OF CODE --------------