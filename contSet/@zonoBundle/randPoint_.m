function p = randPoint_(zB,N,type,varargin)
% randPoint_ - generates a random point within a zonotope bundle
%
% Syntax:
%    p = randPoint_(zB)
%    p = randPoint_(zB,N)
%    p = randPoint_(zB,N,type)
%    p = randPoint_(zB,'all','extreme')
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
% See also: contSet/randPoint, conZonotope/randPoint_

% Authors:       Matthias Althoff
% Written:       18-August-2016 
% Last update:   19-August-2022 (MW, integrate standardized pre-processing)
% Last revision: 27-March-2023 (MW, rename randPoint_)

% ------------------------------ BEGIN CODE -------------------------------

% return all extreme points 
if ischar(N) && strcmp(N,'all')
    p = vertices(zB); return;
end

% generate random points
if strcmp(type,'standard')

    % compute vertices
    Vmat = vertices(zB);
    nrOfVertices = size(Vmat,2);
    if nrOfVertices == 0
        p = Vmat; return
    end

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

        % compute farthest point in this direction 
        [~,x] = supportFunc_(temp,d,'upper');
        p(:,i) = x + c;     
    end
    
else
    throw(CORAerror('CORA:noSpecificAlg',type,zB));
end

% ------------------------------ END OF CODE ------------------------------
