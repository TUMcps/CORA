function P = mptPolytope(pZ,varargin)
% mptPolytope - computes an enclosing polytope for the polynomial zonotope
%
% Syntax:  
%    P = mptPolytope(pZ)
%    P = mptPolytope(pZ,'template',dirs)
%
% Inputs:
%    pZ - polyZonotope object
%    dirs - directions for the template polyhedron
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    pZ = polyZonotope([0;0],[1 0 -1 1;1 1 1 1],[],[1 0 1 2;0 1 1 0]);
%
%    P1 = mptPolytope(pZ);
%    P2 = mptPolytope(pZ,'template',[1 0 1 1;0 1 1 -1]);
%    
%    figure; hold on; xlim([-3,4]); ylim([-2,5]);
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(P1,[1,2],'b');
%    plot(P2,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, zonotope

% Author:       Niklas Kochdumper
% Written:      10-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check if a template polyhedron should be computed
if nargin >= 3 && strcmp(varargin{1},'template') 
    
    % initialize variables
    dir = varargin{2};
    
    C = [dir';dir'];
    d = zeros(size(C,1),1);
    
    counter = 1;
    
    % loop over all directions
    for i = 1:size(dir,2)
        
        dirTemp = dir(:,i)./norm(dir(:,i));
        I = supportFunc_(pZ,dirTemp,'range','split',8,1e-3);
        
        C(counter,:) = dirTemp';
        d(counter) = supremum(I);
        
        C(counter+1,:) = -dirTemp';
        d(counter+1,:) = -infimum(I);
        
        counter = counter + 2;
    end
 
    P = mptPolytope(C,d);
    
else

    Tol = 1e-14;

    % over-approximate all terms where factors with exponents greater than
    % one appear with a zonotope
    ind = find(max(pZ.expMat,[],1) > 1);

    expMat_ = pZ.expMat(:,ind);
    G_ = pZ.G(:,ind);
    c_  = zeros(size(pZ.c));

    pZ.expMat(:,ind) = [];
    pZ.G(:,ind) = [];

    Z = zonotope(polyZonotope(c_,G_,[],expMat_));

    % add over-approximation to the polynomial zonotope
    pZ.c = pZ.c + center(Z);
    pZ.Grest = [pZ.Grest, generators(Z)];

    % determine all potential vertices and remove redundant points
    points = randPoint_(pZ,'all','extreme');
    points = sortrows(points')';

    points_ = zeros(size(points));
    points_(:,1) = points(:,1);
    counter = 1;

    for i = 2:size(points,2)
        if ~all(withinTol(points(:,i),points_(:,counter),Tol))
            counter = counter + 1;
            points_(:,counter) = points(:,i); 
        end
    end

    points = points_(:,1:counter);

    % determine vertices with the n-dimensional convex hull
    ind = convhulln(points');
    ind = unique(squeeze(ind));
    vert = points(:,ind);

    % construct the resulting polytope
    P = mptPolytope(vert');

end

%------------- END OF CODE --------------