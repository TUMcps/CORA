function Z = zonotope(C,varargin)
% zonotope - Over-approximates a capsule by a zonotope
%
% Syntax:  
%    Z = zonotope(C)
%    Z = zonotope(C,order)
%
% Inputs:
%    C - capsule object
%    order - zonotope order of the resulting zonotope
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    C = capsule([0;0],[-2;2],2);
%    Z1 = zonotope(C,2);
%    Z2 = zonotope(C,5);
%
%    figure; hold on
%    plot(C,[1,2],'k');
%    plot(Z1,[1,2],'g');
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Niklas Kochdumper
% Written:      20-Nov-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default input arguments
order = setDefaultValues({5},varargin{:});

% parse input arguments
if nargin >= 2
    order = varargin{1}; 
end

% compute zonotope enclosing the hypersphere
n = length(C.c);
m = order*n-1;

Z = zonotope(ellipsoid(C.r^2*eye(n)),m,'outer:norm');

% constuct enclosing zonotope object
if isempty(C.g)
    Z = Z + C.c; 
else
    Z = zonotope([C.c,C.g,Z.Z(:,2:end)]); 
end

%------------- END OF CODE --------------