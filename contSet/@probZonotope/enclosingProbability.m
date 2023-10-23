function eP = enclosingProbability(probZ,m,dimensions)
% enclosingProbability - Computes the enclosing probability of a 
%    probabilistic zonotope
%
% Syntax:
%    eP = enclosingProbability(probZ,m,dimensions)
%
% Inputs:
%    probZ - probabilistic zonotope object
%    m - m of the mSigma operator
%    dimensions - dimensions on which probability is computed
%
% Outputs:
%    eP - enclosing probability
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    m = 2;
%    dimensions = [1,2];
%    eP = enclosingProbability(probZ,m,dimensions)
%
% Other m-files required: vertices, polytope
% Subfunctions: none

% MAT-files required: none
%
% See also: interval,  vertices

% Authors:       Matthias Althoff
% Written:       08-August-2007
% Last update:   24-August-2007
%                27-August-2007
%                20-March-2015
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

gridpoints=20;

%get Sigma 
origSigma=sigma(probZ);
%project sigma
Sigma(1,1)=origSigma(dimensions(1),dimensions(1));
Sigma(1,2)=origSigma(dimensions(1),dimensions(2));
Sigma(2,1)=origSigma(dimensions(2),dimensions(1));
Sigma(2,2)=origSigma(dimensions(2),dimensions(2));

% fix Sigma for 1D probZ
if all(Sigma(2:4) == 0) && all(probZ.Z(2, :) == 0)
    Sigma(4) = 1;
end
    
%determine mesh size by n-sigma hyperbox
Z = zonotope(probZ,m);
I = interval(Z);
lb = infimum(I);
ub = supremum(I);
x = linspace(lb(dimensions(1)),ub(dimensions(1)),gridpoints);
y = linspace(lb(dimensions(2)),ub(dimensions(2)),gridpoints);

%initialize x-,y- and prob-vector for the mesh
ind=full_fact(1:length(x),1:length(y));
xVector=x(ind(:,1));
yVector=y(ind(:,2)); 
prob=0*xVector;

%check if center of probabilistic zonotope is uncertain
c=center(probZ);
G=probZ.Z(:,2:end);
if isempty(G)
    prob=gaussian([xVector-c(1);yVector-c(2)],Sigma); 
else
    %get uncertain mean
    Z=zonotope(probZ.Z);
    %Compute potential vertices
    V=vertices(Z);
    %Extract projected dimensions
    Vprojected=V(dimensions,:);

    %Plot convex hull of projected vertices
    xPotential=Vprojected(1,:);
    yPotential=Vprojected(2,:);

    %determine vertex indices for convex hull vertices from potential
    %vertices 
    try
        vertexIndices=convhull(xPotential,yPotential);
    
        %Select convex hull vertices
        xCoordinates=xPotential(vertexIndices);
        yCoordinates=yPotential(vertexIndices);
    catch ME
        % ME is thrown if points are collinear (1D)
        [~, max_idx] = max(xPotential);
        [~, min_idx] = min(xPotential);

        xCoordinates=xPotential([min_idx, max_idx]);
        yCoordinates=yPotential([min_idx, max_idx]);
    end

    %points for mean values
    points=[];
    for i=2:length(xCoordinates)
        newX=linspace(xCoordinates(i-1),xCoordinates(i),10);
        newY=linspace(yCoordinates(i-1),yCoordinates(i),10);
        points=[points,[newX;newY]];
    end

    %find inside points
    %get polytope of Z
    P=polytope(Vprojected);
    xyVector = [xVector;yVector];
    insidePoint = xyVector(:,contains(P,xyVector));
    points=[points,insidePoint];

    %for each mean point
    for j=1:size(points,2)
        %get gaussian distribution from the current mean
        c=points(:,j);
        tempP=gaussian([xVector-c(1);yVector-c(2)],Sigma);
        %save maximum values
        prob=max(prob,tempP);
    end
end

%Build probability matrix from probability vector and save x and y values
eP.P=reshape(prob,length(y),length(x));
eP.X=reshape(xVector,length(y),length(x));
eP.Y=reshape(yVector,length(y),length(x));

% ------------------------------ END OF CODE ------------------------------
