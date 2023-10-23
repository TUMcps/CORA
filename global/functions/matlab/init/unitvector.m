function v = unitvector(i,n)
% unitvector - returns the i-th standard unit vector of dimension n
%
% Syntax:
%    v = unitvector(i,n)
%
% Inputs:
%    i - i-th entry is 1
%    n - dimension
%
% Outputs:
%    v - standard unit vector
%
% Example:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Tobias Ladner
% Written:       27-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% omit checks for performance 

if n==0
    % always return empty vector
    v = [];
else
    % init vector of length n
    v = zeros(n,1);
    % set i-th entry to 1
    v(i) = 1;
end

% ------------------------------ END OF CODE ------------------------------
