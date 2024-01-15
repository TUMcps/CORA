function Z = empty(n)
% empty - instantiates an empty zonotope
%
% Syntax:
%    Z = zonotope.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    Z - empty zonotope
%
% Example: 
%    Z = zonotope.empty(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if n <= 0
    throw(CORAerror('CORA:wrongValue','first','positive'));
end

Z = zonotope(zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
