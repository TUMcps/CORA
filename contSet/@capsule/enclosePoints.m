function C = enclosePoints(points)
% enclosePoints - encloses a point cloud by a capsule
%
% Syntax:
%    C = enclosePoints(points)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%
% Outputs:
%    C - capsule object
%
% Example: 
%    p = [-1  1 2 3 2  1 -2 -4 3;...
%          2 -1 3 4 1 -1  2 -3 2];
%    C = capsule.enclosePoints(p);
%    
%    figure; hold on;
%    plot(C);
%    plot(p(1,:),p(2,:),'.k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension of resulting object
[n,nrPoints] = size(points);

% special handling if no or only one point given
if nrPoints == 0
    C = capsule.empty(n);
    return
elseif nrPoints == 1
    C = capsule(points,zeros(n,1),0);
    return
end

% special handling for 1D case
if n == 1
    % center is mean, generator is distance to max value, radius always 0
    c = mean(points);
    g = max(points) - c;
    r = 0;

    % init capsule
    C = capsule(c,g,r);
    return
end


% compute mean of points
c = mean(points,2);

% compute point which is farthest away
dist = vecnorm(points - c);
[~,idxMax] = max(dist,[],2);

% compute generator
g = points(:,idxMax) - c;
g_unit = g ./ vecnorm(g);

% set radius to enclose maximum distance perpendicular to line spanned by
% center +/- generator, we use the formula:
%    distance(a + t*n,p) = || (p-a) - ((p-a)'*n)*n ||_2,
% where the line a+t*n is given by a point 'a', a parameter 't' and the
% vector 'n'; and a point 'p'
r = max(vecnorm(((points-c)' * g_unit)' .* g_unit));

% instantiate capsule
C = capsule(c,g,r);

% ------------------------------ END OF CODE ------------------------------
