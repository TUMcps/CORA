function n = dim(P)
% dim - returns the dimension of the ambient space of a polytope
%
% Syntax:
%    n = dim(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    n - dimension of the polytope
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

if P.isHRep.val
    % either constraints A*x <= b  or  Ae*x == be  given
    if ~isempty(P.A_.val)
        n = size(P.A_.val,2);
    elseif ~isempty(P.Ae_.val)
        n = size(P.Ae_.val,2);
    else
        % constraints, such as zeros(0,n) given
        n = max([size(P.A_.val,2),size(P.Ae_.val,2)]);
    end
elseif P.isVRep.val
    n = size(P.V_.val,1);
end

% ------------------------------ END OF CODE ------------------------------
