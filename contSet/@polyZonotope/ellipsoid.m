function E = ellipsoid(pZ,varargin)
% ellipsoid - computes an enclosing ellipsoid for the polynomial zonotope
%
% Syntax:
%    E = ellipsoid(pZ)
%    E = ellipsoid(pZ,method)
%
% Inputs:
%    pZ - polyZonotope object
%    method - range bounding method that is applied ('interval', ...
%             'split', 'bnb', 'bnbAdv', 'globOpt' or 'bernstein')
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1 1;0 2 1 2],[],[1 0 3 1;0 1 0 2]);
%    E = ellipsoid(pZ);
%
%    figure; hold on;
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(E,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: supportFunc, interval, polytope

% Authors:       Niklas Kochdumper
% Written:       02-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
method = setDefaultValues({'bernstein'},varargin);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {method,'str',{'interval','split','bnb','bnbAdv',...
                    'globOpt','bernstein'}}});

% compute basis with Principal Component Analysis
G = [pZ.G,pZ.GI];
[B,~,~] = svd([-G,G]); 

% compute interval enclosure in the transformed space
I = interval(B'*pZ);
pZ = pZ + (-B*center(I));

% compute Q matrix of the ellipsoid
E = ellipsoid(diag(rad(I).^2));
E = B * E;

% comptue quadratic map of the polynomial zonotope
Q{1} = inv(E.Q);
pZ_ = quadMap(pZ,Q);

% use range bounding to get radius of the ellipsoid
r = supportFunc_(pZ_,1,'upper',method,8,1e-3);

% construct final ellipsoid
E = ellipsoid((E.Q)*r,B*center(I));
    
% ------------------------------ END OF CODE ------------------------------
