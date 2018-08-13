function [P] = mtimes(matrix,P)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
% a pplPolytope
%
% Syntax:  
%    [P] = mtimes(matrix,P)
%
% Inputs:
%    matrix - numerical matrix
%    P - pplPolytope object 
%
% Outputs:
%    P - pplPolytope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      19-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%numeric matrix
if isnumeric(matrix)
    %only C matrix has to be updated
    P.C=P.C*pinv(matrix);
end

%------------- END OF CODE --------------