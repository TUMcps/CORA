function c = center(cZ)
% center - returns a point inside the constrained zonotope. The point is
%    constructed from the chebychev-center of the polytope in the zonotope
%    factor space
%
% Syntax:  
%    c = center(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    c - point inside the constrained zonotope
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZono = conZonotope(Z,A,b);
% 
%    c = center(cZono);
%
%    figure; hold on;
%    plot(cZono,[1,2],'r');
%    plot(c(1),c(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      16-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(cZ.A)
    % no constraints -> zonotope center
    
    if ~isempty(cZ.Z)
        c = cZ.Z(:,1);
    else
        c = [];
    end

else
    % constraints -> compute chebychev center
        
    % construct inequality constraints for the unit cube
    n = size(cZ.Z,2)-1;
    A = [eye(n);-eye(n)];
    b = [ones(n,1);ones(n,1)];
    
    % calculate null space of the constraints
    Neq = null(cZ.A); 

    % if null space is empty, there is only one solution for factors
    % -> constrained zonotope is a single point
    if isempty(Neq)
        c = cZ.Z(:,1) + sum(cZ.Z(:,2:end) * diag(cZ.A \ cZ.b),2);
        return
    end
    
    % Calculate a single point that satisfies the constraints
    x0 = pinv(cZ.A)*cZ.b;
    
    % transform the constraints to the null space
    A_ = A*Neq;
    b_ = b-A*x0;
    
    % create mptPolytope
    poly = mptPolytope(A_,b_);
    
    % compute chebychev center in the zonotope-factor null-space
    P = get(poly,'P');
    c = chebyCenter(P);
    
    % convert center back to the normal zonotope factor space
    c_ = Neq*c.x + x0;
    
    % compute center of the constraint zonotope using the the factors
    % from the chebychev center in the factor space
    c = cZ.Z(:,1) + cZ.Z(:,2:end) * c_;
end

%------------- END OF CODE --------------