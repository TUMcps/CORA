function res = zonoBundle(obj)
% conZonotope - convert an conZonotope object into a zonotope bundle object
%
% Syntax:  
%    res = zonoBundle(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%    res - zonoBundle object
%
% Example:
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZ = conZonotope(Z,A,b);
%
%    zB = zonoBundle(cZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope/zonoBundle

% Author:       Niklas Kochdumper
% Written:      26-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    if isempty(obj.A)
        
        res = zonoBundle(zonotope(obj.Z)); 
       
    else
        
        % compute point satisfying all constraints with pseudo inverse
        p_ = pinv(obj.A)*obj.b;

        % compute null-space of constraints
        T = null(obj.A);

        % transform boundary constraints of the factor hypercube
        m = size(obj.A,2);

        A = [eye(m);-eye(m)];
        b = ones(2*m,1);

        A_ = A*T;
        b_ = b - A*p_;

        % enclose null-space polytope by a zonotope bundle
        poly_ = mptPolytope(A_,b_);
        zB_ = zonoBundle(poly_);

        % transform back to original space
        c = obj.Z(:,1);
        G = obj.Z(:,2:end);
        
        res = c + G * p_ + (G * T) * zB_;     
        
    end

%------------- END OF CODE --------------