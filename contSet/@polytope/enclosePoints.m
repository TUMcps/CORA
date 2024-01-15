function P_out = enclosePoints(points)
% enclosePoints - encloses a point cloud with a polytope
%
% Syntax:
%    P_out = enclosePoints(points)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%
% Outputs:
%    P_out - polytope object
%
% Example: 
%    mu = [2 3];
%    sigma = [1 1.5; 1.5 3];
%    points = mvnrnd(mu,sigma,100)';
%
%    P = polytope.enclosePoints(points);
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(P,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclosePoints, interval/enclosePoints

% Authors:       Niklas Kochdumper, Victor Gassmann
% Written:       05-May-2020
% Last update:   17-March-2021 (now also works for degenerate point cloud)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

r = rank(points);
[n,N] = size(points);
Q = eye(n);

% check whether points are lower-dimensional
if rank(points)<size(points,1)
    [Q,~] = qr(points);
    points = Q'*points;
    % remove zeros
    points(r+1:end,:) = [];
end

% compute convex hull
if size(points,1)>1
    ind = convhulln(points');
else
    [~,ii_max] = max(points);
    [~,ii_min] = min(points);
    ind = [ii_max,ii_min];
end

ind = reshape(ind,[numel(ind),1]);
ind = unique(ind);

points = points(:,ind);

% add zeros again and backtransform
points(end+1:end+n-r,:) = zeros(n-r,length(ind));
points = Q*points;

% construct polytope object (H-representation computed in constructor)
P_out = polytope(points);

% ------------------------------ END OF CODE ------------------------------
