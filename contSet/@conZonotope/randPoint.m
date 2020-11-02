function p = randPoint(cZ)
% randPoint - generates a random point inside a constrained zonotope
%
% Syntax:  
%    p = randPoint(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    cZ = conZonotope.generateRandom(2);
%
%    points = zeros(2,100);
%    for i = 1:100
%       points(:,i) = randPoint(cZ);
%    end
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPointExtreme

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(cZ.A)   

           % no constraints -> call superclass method
           p = randPoint(zonotope(cZ.Z));

    else                

           % construct inequality constraints for the unit cube
           n = size(cZ.Z,2)-1;
           A = [eye(n);-eye(n)];
           b = [ones(n,1);ones(n,1)];

           % calculate null space of the constraints
           Neq = null(cZ.A); 

           % calculate a single point that satisfies the constraints
           x0 = pinv(cZ.A)*cZ.b;

           % transform the constraints to the null space
           A_ = A*Neq;
           b_ = b-A*x0;

           % create mptPolytope
           poly = mptPolytope(A_,b_);

           % compute chebychev center in the zonotope-factor null-space
           p = randPoint(poly);

           % convert center back to the normal zonotope factor space
           p_ = Neq*p + x0;

           % compute center of the constraint zonotope using the the
           % factors from the chebychev center in the factor space
           p = cZ.Z(:,1) + cZ.Z(:,2:end) * p_;
    end
end

%------------- END OF CODE ----------