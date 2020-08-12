function poly = enclosePoints(points)
% enclosePoints - enclose a point cloud with a polytope
%
% Syntax:  
%    poly = enclosePoints(points)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%
% Outputs:
%    poly - mptPolytope object
%
% Example: 
%    mu = [2 3];
%    sigma = [1 1.5; 1.5 3];
%    points = mvnrnd(mu,sigma,100)';
%
%    poly = mptPolytope.enclosePoints(points);
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(poly,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclosePoints, interval/enclosePoints

% Author: Niklas Kochdumper
% Written: 05-May-2020
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % compute convex hull
    ind = convhulln(points');
    
    ind = reshape(ind,[numel(ind),1]);
    ind = unique(ind);
    
    points = points(:,ind);
    
    % construct mptPolytope object
    poly = mptPolytope(points');

%------------- END OF CODE --------------