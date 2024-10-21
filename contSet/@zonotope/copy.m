function Z_out = copy(Z)
% copy - copies the zonotope object (used for dynamic dispatch)
%
% Syntax:
%    Z_out = copy(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z_out - copied zonotope object
%
% Example: 
%    Z = zonotope([1;0],[1 0 -1; 0 1 1]);
%    Z_out = copy(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
Z_out = zonotope(Z);

% ------------------------------ END OF CODE ------------------------------
