function p = randPoint(zB,varargin)
% randPoint - generates a random point within a zonotope bundle
%
% Syntax:  
%    p = randPoint(zB)
%    p = randPoint(zB,N)
%    p = randPoint(zB,N,type)
%    p = randPoint(zB,'all','extreme')
%
% Inputs:
%    zB - zonoBundle object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
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
% Last update:  19-August-2022 (MW, integrate standardized pre-processing)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_randPoint('zonoBundle',zB,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of vars
    p = vars{1}; return
else
    % assign variables
    zB = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
    % sampling with 'gaussian' is done in contSet method
    if strcmp(type,'gaussian')
        p = randPoint@contSet(zB,N,type,pr); return
    end
end


% return all extreme points 
if ischar(N) && strcmp(N,'all')
    p = vertices(zB); return;
end

% generate random points
if strcmp(type,'standard')

    % compute vertices
    Vmat = vertices(zB);
    nrOfVertices = length(Vmat(1,:));

    % random convex combination
    alpha = rand(nrOfVertices,N);
    alphaNorm = alpha./sum(alpha,1);

    % random points
    p = Vmat * alphaNorm;
    
elseif strcmp(type,'extreme')
    
    p = zeros(dim(zB),N);
    
    % center polytope at origin
    c = center(zB);
    temp = zB + (-c);
    
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
    
end

%------------- END OF CODE --------------