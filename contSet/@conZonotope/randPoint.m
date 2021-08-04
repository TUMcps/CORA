function p = randPoint(obj,varargin)
% randPoint - generates a random point inside a constrained zonotope
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,'all','extreme')
%
% Inputs:
%    obj - conZonotope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    cZ = conZonotope.generateRandom(2);
%
%    points = randPoint(cZ,100);
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'standard';
    if nargin > 1 && ~isempty(varargin{1})
       N = varargin{1}; 
    end
    if nargin > 2 && ~isempty(varargin{2})
       type = varargin{2}; 
    end
    
    % call zonotope method if no constraints are present
    if isempty(obj.A)
       p = randPoint(zonotope(obj.Z),N,type);
       return;
    end
    
    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        p = vertices(obj); return;
    end
    
    % generate random points
    p = zeros(dim(obj),N);
    
    if strcmp(type,'standard')
        for i = 1:N
            p(:,i) = randPointNormal(obj);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            p(:,i) = randPointExtreme(obj);
        end
    else
        [msg,id] = errWrongInput('type');
        error(id,msg);
    end
end


% Auxiliary Functions -----------------------------------------------------

function p = randPointNormal(obj)    
% generate random point within the constrained zonotope

    % construct inequality constraints for the unit cube
    n = size(obj.Z,2)-1;
    A = [eye(n);-eye(n)];
    b = [ones(n,1);ones(n,1)];

    % calculate null space of the constraints
    Neq = null(obj.A); 

    % calculate a single point that satisfies the constraints
    x0 = pinv(obj.A)*obj.b;

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
    p = obj.Z(:,1) + obj.Z(:,2:end) * p_;
end

function p = randPointExtreme(obj)
% generate random point on boundary of a constrained zonotope

    % center constrained zonotope at origin
    c = center(obj);
    obj = obj + (-c);

    % select random direction
    n = length(c);
    d = rand(n,1) - 0.5*ones(n,1);
    d = d./norm(d);

    % compute farest point in this direction that is still in set
    [~,x] = supportFunc(obj,d);
    p = x + c;
end

%------------- END OF CODE ----------