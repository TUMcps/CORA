function val = dim(P)
% dim - returns the dimension of the ambient space of a polytope
%
% Syntax:
%    val = dim(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    val - dimension of the polytope
%
% Example:
%    A = [1 2; -1 2; -2 -2; 1 -2];
%    b = ones(4,1);
%    P = polytope(A,b);
%    dim(P)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       28-March-2022
% Last update:   01-December-2022
%                19-July-2023 (MW, support zeros(0,n) constraints)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% either constraints A*x <= b  or  A*x <= be  given
if ~isempty(P.A)
    val = size(P.A,2);
elseif ~isempty(P.Ae)
    val = size(P.Ae,2);
else
    % constraints, such as zeros(0,n) given
    val = max([size(P.A,2),size(P.Ae,2)]);
end

% ------------------------------ END OF CODE ------------------------------
