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

        % normalize halfspace directions to length 1
        temp = diag(1./sqrt(sum(A_.^2,2)));
        A_ = temp * A_; b_ = temp * b_;
        
        % enclose null-space polytope by a zonotope bundle
        d = sqrt(m);
        p = size(T,2);
        Z = cell(m,1);
        
        for i = 1:m
            
            % compute basis orthogonal to current normal vector
            B = gramSchmidt(A_(i,:)');
            
            % compute interval enclousre in the transformed space
            int = interval([-b_(m+i);-d*ones(p-1,1)],[b_(i);d*ones(p-1,1)]);
            Z{i} = B * zonotope(int);
        end
        
        zB_ = zonoBundle(Z);

        % transform back to original space
        c = obj.Z(:,1);
        G = obj.Z(:,2:end);
        
        res = c + G * p_ + (G * T) * zB_;     
        
    end

%------------- END OF CODE --------------