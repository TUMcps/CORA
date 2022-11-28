function res = isequal(C1,C2,varargin)
% isequal - checks if two capsules are equal
%
% Syntax:  
%    res = isequal(C1,C2)
%
% Inputs:
%    C1 - capsule object
%    C2 - capsule object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
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
%               01-June-2022 (MW, correct dimension mismatch case)
% Last revision:---

%------------- BEGIN CODE --------------

% default values
tol = setDefaultValues({eps},varargin{:});

% check input arguments
inputArgsCheck({{C1,'att','capsule'};
                {C2,'att','capsule'};
                {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

if dim(C1) ~= dim(C2)
    res = false; return
end

% check for equality
res = all(abs(center(C1) - center(C2)) < tol) && ... % center
    all(abs(C1.g - C2.g) < tol) && ... % generator
    abs(C1.r - C2.r) < tol; % radius

%------------- END OF CODE --------------