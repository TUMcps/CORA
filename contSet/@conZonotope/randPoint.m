function p = randPoint(cZ,varargin)
% randPoint - generates a random point inside a constrained zonotope
%
% Syntax:  
%    p = randPoint(cZ)
%    p = randPoint(cZ,N)
%    p = randPoint(cZ,N,type)
%    p = randPoint(cZ,'all','extreme')
%
% Inputs:
%    cZ - conZonotope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    cZ = conZonotope.generateRandom('Dimension',2);
%    p = randPoint(cZ,100);
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(p(1,:),p(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  19-August-2022 (MW, integrate standardized pre-processing)
% Last revision:---

%------------- BEGIN CODE --------------

    % pre-processing
    [res,vars] = pre_randPoint('conZonotope',cZ,varargin{:});
    
    % check premature exit
    if res
        % if result has been found, it is stored in the first entry of vars
        p = vars{1}; return
    else
        % assign variables
        cZ = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
        % sampling with 'gaussian' is done in contSet method
        if strcmp(type,'gaussian')
            p = randPoint@contSet(cZ,N,type,pr); return
        end
    end

    % call zonotope method if no constraints are present
    if isempty(cZ.A)
        p = randPoint(zonotope(cZ.Z),N,type,pr); return
    end
    
    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        p = vertices(cZ); return
    end
    
    % generate random points
    p = zeros(dim(cZ),N);
    
    if strcmp(type,'standard')
        for i = 1:N
            p(:,i) = randPointNormal(cZ);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            p(:,i) = randPointExtreme(cZ);
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function p = randPointNormal(cZ)    
% generate random point within the constrained zonotope

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

    % instantiate mptPolytope
    poly = mptPolytope(A_,b_);

    % compute Chebychev center in the zonotope-factor null-space
    p = randPoint(poly);

    % convert center back to the normal zonotope factor space
    p_ = Neq*p + x0;

    % compute center of the constraint zonotope using the factors from the
    % Chebychev center in the factor space
    p = cZ.Z(:,1) + cZ.Z(:,2:end) * p_;
end

function p = randPointExtreme(cZ)
% generate random point on boundary of a constrained zonotope

    % center constrained zonotope at origin
    c = center(cZ);
    cZ = cZ + (-c);

    % select random direction
    n = dim(cZ);
    d = rand(n,1) - 0.5*ones(n,1);
    d = d./norm(d);

    % compute farest point in this direction that is still in set
    [~,x] = supportFunc(cZ,d);
    p = x + c;
end

%------------- END OF CODE ----------