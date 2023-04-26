function r = is_close(A, B)
% is_close - test if two matrices are identical up to a tolerance
%
%
% Syntax:
%    res = is_close(A, B)
%
% Inputs:
%    A - matrix
%    B - matrix (same dimension)
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Tobias Ladner
% Written:      24-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

r = max(max(abs(A-B))) <= 1.0e-12;
m = max(max(abs([A, B])));
if m == 0
    m = 1;
end
r = max(max(abs(A-B))) / m <= 1.0e-10;
end

%------------- END OF CODE --------------