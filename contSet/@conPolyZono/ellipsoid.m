function E = ellipsoid(cPZ,varargin)
% ellipsoid - computes an enclosing ellipsoid of a constrained polynomial 
%    zonotope
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
%    E = [1 0 2 0; 0 1 1 0; 0 0 0 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [1 0 0; 0 1 2; 0 1 0];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    E = ellipsoid(cPZ);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(E,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: supportFunc, interval, zonotope

% Authors:       Niklas Kochdumper
% Written:       05-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
method = setDefaultValues({'split'},varargin); 

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {method,'str',{'interval','split','quadProg',...
                    'conZonotope','interval'}}});

% compute basis with Principal Component Analysis
G = [cPZ.G,cPZ.GI];
[B,~,~] = svd([-G,G]); 

% compute interval enclosure in the transformed space
I = interval(B'*cPZ);
cPZ = cPZ + (-B*center(I));

% compute Q matrix of the ellipsoid
E = ellipsoid(diag(rad(I).^2));
E = B * E;

% comptue quadratic map of the polynomial zonotope
Q{1} = inv(E.Q);
cPZ_ = quadMap(cPZ,Q);

% use range bounding to get radius of the ellipsoid
r = supportFunc_(cPZ_,1,'upper',method,8);

% construct final ellipsoid
E = ellipsoid((E.Q)*r,B*center(I));
    
% ------------------------------ END OF CODE ------------------------------
