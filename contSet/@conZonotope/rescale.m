function res = rescale(obj,varargin)
% resclae - Rescales the domains for the factors ksi of a constrained
%           zonotope
%
% Syntax:  
%    res = rescale(obj)
%    res = rescale(obj, method)
%
% Inputs:
%    obj - c-zonotope object
%    method - method used to determine the tighend domain for the zonotope
%             factors ('exact' or 'iter')
%
% Outputs:
%    res - c-zonotope object
%
% Example:
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    cRes = rescale(cZono,'iter');
%    
%    hold on
%    plotZono(cZono,[1,2]);
%    plotZono(cRes,[1,2],{'g'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      18-December-2017
% Last update:  12-May-2018
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
method = 'exact';
if nargin == 2
   method = varargin{1};
end

% determine the new domain for each factor ksi
if isempty(obj.A)
    error('No constrains left')
elseif isempty(obj.ksi)
    if strcmp(method,'iter')
        [domKsi, R] = ksi_iterative( obj );
    else
        [domKsi,~] = ksi_optimizer( obj );
    end
    
    ksi_l = infimum(domKsi);
    ksi_u = supremum(domKsi);
else
    ksi_l = min(obj.ksi,[],2);
    ksi_u = max(obj.ksi,[],2);
end

% new basis
ksi_m = (ksi_u + ksi_l)/2;
ksi_r = (ksi_u - ksi_l)/2;

% rescale R (Equation (A.3) in reference paper [1])
if strcmp(method,'iter')
    rho_l = infimum(R);
    rho_u = supremum(R);

    obj.R = [(rho_l - ksi_m)./ksi_r, (rho_u - ksi_m)./ksi_r];
end

% rescale c-zonotope (Equation (24) in reference paper [1])
temp = diag(ksi_r);

G = obj.Z(:, 2:end);
c = obj.Z(:, 1) + G*ksi_m;
obj.Z = [c, G * temp];
obj.b = obj.b - obj.A * ksi_m;
obj.A = obj.A * temp;

if ~isempty(obj.ksi)
   temp = ones(1,size(obj.ksi,2));
   ksi = obj.ksi - ksi_m*temp;
   obj.ksi = zeros(size(ksi)) + ksi.*((1./ksi_r)*temp);
end

% output the result
res = obj;

end