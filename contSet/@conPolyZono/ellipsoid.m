function E = ellipsoid(cPZ,varargin)
% ellipsoid - computes an enclosing ellipsoid of a constrained polynomial 
%             zonotope
%
% Syntax:  
%    E = ellipsoid(cPZ)
%    E = ellipsoid(cPZ,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - range bounding method that is applied ('interval', ...
%             'split', 'quadProg', 'conZonotope', or 'interval')
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    c = [2;-1];
%    G = [2 1 2 1; 0 2 2 1];
%    expMat = [1 0 2 0; 0 1 1 0; 0 0 0 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    expMat_ = [1 0 0; 0 1 2; 0 1 0];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    E = ellipsoid(cPZ);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%    plot(E,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: supportFunc, interval, zonotope

% Author:       Niklas Kochdumper
% Written:      05-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    method = 'split';
    if nargin > 1
       method = varargin{1}; 
    end

    % compute basis with Principal Component Analysis
    G = [cPZ.G,cPZ.Grest];
    [B,~,~] = svd([-G,G]); 
    
    % compute interval enclosure in the transformed space
    int = interval(B'*cPZ);
    cPZ = cPZ + (-B*center(int));
    
    % compute Q matrix of the ellipsoid
    E = ellipsoid(diag(rad(int).^2));
    E = B * E;
    
    % comptue quadratic map of the polynomial zonotope
    Q{1} = inv(E.Q);
    cPZ_ = quadMap(cPZ,Q);
    
    % use range bounding to get radius of the ellipsoid
    r = supportFunc(cPZ_,1,'upper',method);
    
    % construct final ellipsoid
    E = ellipsoid((E.Q)*r,B*center(int));
    
end
    
%------------- END OF CODE --------------