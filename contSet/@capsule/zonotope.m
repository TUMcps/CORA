function Z = zonotope(obj,varargin)
% zonotope - Over-approximate a capsule by a zonotope
%
% Syntax:  
%    Z = zonotope(obj)
%    Z = zonotope(obj,order)
%
% Inputs:
%    obj - capsule object
%    order - zonotope order of the resulting zonotope
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    C = capsule([0;0],[-2;2],2);
%
%    zono1 = zonotope(C,2);
%    zono2 = zonotope(C,5);
%
%    figure
%    hold on
%    plot(C,[1,2],'k');
%    plot(zono1,[1,2],'g');
%    plot(zono2,[1,2],'r');
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

    % parse input arguments
    order = 5;
    
    if nargin >= 2 && ~isempty(varargin{1})
       order = varargin{1}; 
    end
    
    % compute zonotope enclosing the hypersphere
    n = length(obj.c);
    m = order*n-1;
    
    zono = zonotope(ellipsoid(obj.r^2*eye(n)),m,'o:norm');
    
    % constuct enclosing zonotope object
    if isempty(obj.g)
       Z = zono + obj.c; 
    else
       Z = zonotope([obj.c,obj.g,zono.Z(:,2:end)]); 
    end

%------------- END OF CODE --------------