function p = randPoint(obj,varargin)
% randPoint - generates a random point within a zonotope bundle
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,'all','extreme')
%
% Inputs:
%    obj - zonotope bundle object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    points = randPoint(zB,100);
%   
%    figure; hold on;
%    plot(zB,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/randPoint

% Author:       Matthias Althoff
% Written:      18-August-2016 
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

    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        p = vertices(obj); return;
    end

    % generate random points
    if strcmp(type,'standard')

        % compute vertices
        Vmat = vertices(obj);
        nrOfVertices = length(Vmat(1,:));

        % random convex combination
        alpha = rand(nrOfVertices,N);
        alphaNorm = alpha./sum(alpha,1);

        % random points
        p = Vmat * alphaNorm;
        
    elseif strcmp(type,'extreme')
        
        p = zeros(dim(obj),N);
        
        % center polytope at origin
        c = center(obj);
        temp = obj + (-c);
        
        % loop over all points
        for i = 1:N

            % select random direction
            n = length(c);
            d = rand(n,1) - 0.5*ones(n,1);
            d = d./norm(d);

            % compute farest point in this direction 
            [~,x] = supportFunc(temp,d);
            p(:,i) = x + c;     
        end
        
    else
        [msg,id] = errWrongInput('type');
        error(id,msg);
    end
end

%------------- END OF CODE --------------