function zB = zonoBundle(cZ)
% zonoBundle - converts an constrained zonotope to a zonotope bundle
%
% Syntax:
%    zB = zonoBundle(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    zB - zonoBundle object
%
% Example:
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%
%    zB = zonoBundle(cZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/zonoBundle

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
if isempty(cZ.A)
    
    zB = zonoBundle(zonotope(cZ.c,cZ.G)); 
   
else
    
    % compute point satisfying all constraints with pseudo inverse
    p_ = pinv(cZ.A)*cZ.b;

    % compute null-space of constraints
    T = null(cZ.A);

    % transform boundary constraints of the factor hypercube
    m = size(cZ.A,2);

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
    c = cZ.c;
    G = cZ.G;
    
    zB = c + G * p_ + (G * T) * zB_;     
    
end

% ------------------------------ END OF CODE ------------------------------
