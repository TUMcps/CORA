function R = uminus(R)
% uminus - Overloaded the unary '-' operator
%
% Syntax:
%    R = -R
%    R = uminus(R)
%
% Inputs:
%    R - reachSet object 
%
% Outputs:
%    R - transformed reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/uminus

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

R = -1 * R;

% ------------------------------ END OF CODE ------------------------------
