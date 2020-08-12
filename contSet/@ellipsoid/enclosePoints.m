function E = enclosePoints(points)
% enclosePoints - enclose a point cloud with an ellipsoid
%
% Syntax:  
%    E = enclosePoints(points)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    mu = [2 3];
%    sigma = [1 1.5; 1.5 3];
%    points = mvnrnd(mu,sigma,100)';
%
%    E = ellipsoid.enclosePoints(points);
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(E,[1,2],'r');
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

    % compute the arithmetic mean of the points
    mean = sum(points,2)/length(points(1,:));

    % obtain sampling matrix
    translation = mean*ones(1,length(points(1,:)));
    sampleMatrix = points-translation;

    % compute the covariance matrix
    C = cov(sampleMatrix');

    % singular value decomposition
    [U,~,~] = svd(C);
    
    % required ellipsoid radius in transformed space
    orientedMatrix=U'*sampleMatrix;
    m1 = max(orientedMatrix,[],2);
    m2 = min(orientedMatrix,[],2);

    r = max([m1,m2],[],2);
    
    % enclosing ellipsoid in the transformed space
    Q = diag(r.^2);
    
    maxDist = max(sum(diag(1./r.^2)*orientedMatrix.^2,1));
    Q = Q * maxDist;
    
    E = ellipsoid(Q);
    
    % transform back to original space
    E = U*E  + mean;

%------------- END OF CODE --------------