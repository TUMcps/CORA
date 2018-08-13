function V = potVertices(obj)
% potVertices - calculate all potential vertices of a constrained zonotope
%               object. The points vertices are either real vertices or are
%               located inside the constrained zonotope 
%
% Syntax:  
%    res = potVertices(obj)
%
% Inputs:
%    obj - matrix of size (n,m) containing the m potential vertices
%
% Outputs:
%    res - c-zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% Calculate extreme points for the zonotope factors ksi
if isempty(obj.ksi)

    % remove all constraints for which the elimination does not result
    % in an over-approximation (less ksi-dimensions -> speed-up)
    obj = reduce(obj,'redConstr');
    
    % bounding constraints
    n = size(obj.Z,2)-1;
    A = [diag(ones(n,1));-diag(ones(n,1))];
    b = ones(2*n,1);
    
    % calculate the extreme points in ksi-space (= polytope vertices)
    ksi = lcon2vert(A,b,obj.A,obj.b);
    obj.ksi = ksi';

end
    
% Calculate the corresponding zonotope points (ksi-space -> real space)
V = obj.Z * [ones(1,size(obj.ksi,2)); obj.ksi];

%------------- END OF CODE --------------