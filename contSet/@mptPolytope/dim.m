function n = dim(P)
% dim - returns the dimension of the ambient space of a polytope
%
% Syntax:  
%    n = dim(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      16-July-2015
% Last update:  03-February-2017
%               17-March-2017
% Last revision:---

%------------- BEGIN CODE --------------

%return dimension
try %MPT3
    n = P.P.Dim;
catch
    [H,~] = double(P.P);
    n = length(H(end,:));
end

%------------- END OF CODE --------------