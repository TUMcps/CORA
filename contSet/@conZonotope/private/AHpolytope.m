function [G,c,A,b] = AHpolytope(cZ)
% AHpolytope - converts a constrained zonotope to a AH-polytope 
%              x = c + G z with A z <= b
%
% Syntax:  
%    [G,c,A,b] = AHpolytope(cZ)
%
% Inputs:
%    cZ - conZonotope objects
%
% Outputs:
%    c - center vector of AH-polytope
%    G - generator matrix of AH-polytope
%    A - constraint matrix of AH-polytope
%    b - constraint vector of AH-polytope
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or, zonotope/or

% Author:       Niklas Kochdumper
% Written:      14-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % check if set is a zonotope
    if isempty(cZ.A)
        
        c = cZ.Z(:,1);
        G = cZ.Z(:,2:end);
        
        m = size(G,2);
        A= [eye(m);-eye(m)];
        b = ones(2*m,1);
        
    else

        % compute point satisfying all constraints with pseudo inverse
        p_ = pinv(cZ.A)*cZ.b;

        % compute null-space of constraints
        T = null(cZ.A);

        % transform boundary constraints of the factor hypercube
        m = size(cZ.A,2);
        m_ = size(T,2);

        A_ = [eye(m);-eye(m)];
        b_ = ones(2*m,1);

        A = A_*T;
        b = b_ - A_*p_;

        % compute updated center and generator matrix
        c = cZ.Z(:,1);
        G = cZ.Z(:,2:end);

        c = c + G*p_;
        G = G*T;
    end
end

%------------- END OF CODE --------------