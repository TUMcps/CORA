function E = ellipsoid(pZ)
% ellipsoid - computes an enclosing ellipsoid for the polynomial zonotope
%
% Syntax:  
%    E = ellipsoid(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1 1;0 2 1 2],[],[1 0 3 1;0 1 0 2]);
%    E = ellipsoid(pZ);
%
%    figure
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(E,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, mptPolytope

% Author:       Niklas Kochdumper
% Written:      02-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % center polynomial zonotope at origin
    c = center(pZ);
    pZ = pZ + (-c);

    % compute basis with Principal Component Analysis
    G = [pZ.G,pZ.Grest];
    [B,~,~] = svd([-G,G]); 
    
    % compute interval enclosure in the transformed space
    int = interval(B'*pZ);
    
    % compute Q matrix of the ellipsoid
    E = ellipsoid(diag(rad(int).^2));
    E = B * E;
    
    % comptue quadratic map of the polynomial zonotope
    Q{1} = inv(E.Q);
    pZ_ = quadMap(pZ,Q);
    
    % use range bounding to get radius of the ellipsoid
    r = supportFunc(pZ_,1,'upper','bernstein');
    
    % construct final ellipsoid
    E = ellipsoid((E.Q)*r,c);
    
%------------- END OF CODE --------------