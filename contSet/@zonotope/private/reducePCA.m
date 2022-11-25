function Zred = reducePCA(Z,order)
% reducePCA - apply principal component analysis
%
% Syntax:  
%    [Zred,t]=reducePCA(Z)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      18-October-2013
% Last update:  11-October-2017
% Last revision:---

%------------- BEGIN CODE --------------

% initialize Z_red
Zred=Z;

% pick generators to reduce
[center, Gunred, Gred] = pickedGenerators(Z,order);

if ~isempty(Gred)
    %obtain matrix of points from generator matrix
    V = [Gred,-Gred];

    %compute the arithmetic mean of the vertices
    mean=sum(V,2)/length(V(1,:));

    %obtain sampling matrix
    translation=mean*ones(1,length(V(1,:)));
    sampleMatrix=V-translation;

    %compute the covariance matrix
    C=cov(sampleMatrix');

    %singular value decomposition
    [U,~,~] = svd(C);

    % map generators
    Gtrans = U.'*Gred;

    % box generators
    Gbox=diag(sum(abs(Gtrans),2));

    % transform generators back
    Gred = U*Gbox;
end

%build reduced zonotope
Zred.Z=[center,Gunred,Gred];

%------------- END OF CODE --------------
