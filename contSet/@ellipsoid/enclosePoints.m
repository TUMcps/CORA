function E = enclosePoints(points,method)
% enclosePoints - enclose a point cloud with an ellipsoid
%
% Syntax:  
%    E = enclosePoints(points)
%    E = enclosePoints(points,method);
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%    method - method to compute the enclosing ellipsoid
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    mu = [2 3];
%    sigma = [1 1.5; 1.5 3];
%    points = mvnrnd(mu,sigma,100)';
%
%    E1 = ellipsoid.enclosePoints(points);
%    E2 = ellipsoid.enclosePoints(points,'min-vol');
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(E1,[1,2],'r');
%    plot(E2,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclosePoints, interval/enclosePoints

% Author: Niklas Kochdumper, Victor Gassmann
% Written: 05-May-2020
% Last update: 15-March-2021
% Last revision: ---

%------------- BEGIN CODE --------------
if ~exist('method','var')
    method = 'cov';
end

% remove bias
c = mean(points,2);
points = points-c;
n = size(points,1);

% handle degenerate case
[T,S,~] = svd(points);
n_nd = rank(ellipsoid(S(:,1:n)));
points = T'*points;
n_d = n-n_nd;
r_d = zeros(n_d,1);
if n_nd<n
    if n_nd==0
        % all zero matrix
        E = ellipsoid(zeros(n),zeros(n,1));
        return;
    end
    % essentially, this is pca; since we know that the data is
    % "degenerate", the last row has to be ~zero, as the variance is ~0
    % remove zeros 
    % save rest to add later (don't just ignore the very small parts)
    r_d = 1/2*(max(points(n_nd+1:end,:),[],2)-min(points(n_nd+1:end,:),[],2));
    points(n_nd+1:end,:) = [];
end

if strcmp(method,'cov')

    % compute the covariance matrix
    C = cov(points');

    % singular value decomposition
    [U,~,~] = svd(C);
    
    % required ellipsoid radius in transformed space
    orientedMatrix=U'*points;
    m1 = max(orientedMatrix,[],2);
    m2 = min(orientedMatrix,[],2);

    nt = max([m1,m2],[],2);
    
    % enclosing ellipsoid in the transformed space
    Q = diag(nt.^2);
    
    maxDist = max(sum(diag(1./nt.^2)*orientedMatrix.^2,1));
    Q = Q * maxDist;
    
    E = ellipsoid(Q);
    
    % transform back to original space
    E = U*E;
    
elseif strcmp(method,'min-vol')
    nt = size(points,1);
    Q = sdpvar(nt);
    X = points;
    F = cone([ones(size(X,2),1),X'*Q]');
    options = sdpsettings;
    options.verbose = 0;
    % if sdpt3 is installed, use it
    if exist('sdpt3','file')
        options.solver = 'sdpt3';
    end
    % solve optimization problem
    diagnostics = optimize(F, -logdet(Q), options);
    if diagnostics.problem ~= 0
        throw(errOptNumIssue());
    end
    E = ellipsoid(inv(value(Q)^2));
else
    error('Method not supported!');
end

% backtransform and add back center
E_ext = E;
for i=1:n_d
    E_ext = cartProd(E_ext,ellipsoid(r_d(i)^2));
end
E = T*E_ext + c;
%------------- END OF CODE --------------