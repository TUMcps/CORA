function res = isFullDim(obj)
% isFullDim - check if a polytope is full-dimensional
%
% Syntax:  
%    res = isFullDim(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    res - 1 if mptPolytope is full-dimensional, 0 else
%
% Example:
%    poly1 = mptPolytope([1 0;0 1;-1 0; 0 -1],[1;1;0;0]);
%    poly2 = mptPolytope([1 0;0 1;-1 0; 0 -1],[1;0;0;0]);
%
%    isFullDim(poly1)
%    isFullDim(poly2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;

    % get object properties
    A = obj.P.A;
    b = obj.P.b;

    % normalize constraint matrix
    n = sqrt(sum(A.^2,2));
    A = A./n;
    b = b./n;
    
    % loop over all constraints
    for i = 1:length(b)
        
       [res_,ind] = ismembertol(-A(i,:),A,eps,'ByRows',true);
        
       if res_ && b(ind) == -b(i)
           res = 0;
           return;
       end
    end

%------------- END OF CODE --------------