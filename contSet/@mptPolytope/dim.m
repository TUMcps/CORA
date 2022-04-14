function dimension = dim(obj)
% dim - returns the dimension of the mptPolytope
%
% Syntax:  
%    dimension = dim(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    dimension - dimension of the mptPolytope
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
    dimension = length(obj.P.A(end,:));
catch
    [H,~]=double(obj.P);
    dimension = length(H(end,:));
end

%------------- END OF CODE --------------