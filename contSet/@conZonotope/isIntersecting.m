function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if constrained zonotope obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - conZonotope object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    % generate constrained zonotopes
%    Z = [0 2 -2 1;0 1.5 1 -1.5];
%    A = [1 1 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);
% 
%    Z = [1 2 0 0;1 1 1 0];
%    A = [1 1 -1];
%    b = 0;
%    cZono2 = conZonotope(Z,A,b);
%
%    Z = [3 2 0 0;4 1 1 0];
%    A = [1 1 -1];
%    b = 0;
%    cZono3 = conZonotope(Z,A,b);
%
%    % check for intersection
%    isIntersecting(cZono1,cZono2)
%    isIntersecting(cZono1,cZono3)
%
%    % visualization
%    figure; hold on
%    plot(cZono1,[1,2],'b');
%    plot(cZono2,[1,2],'g');
%
%    figure; hold on
%    plot(cZono1,[1,2],'b');
%    plot(cZono3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    % get conZonotope object
    if ~isa(obj1,'conZonotope')
       temp = obj1;
       obj1 = obj2;
       obj2 = temp;
    end
    
    % intersection with halfspace, hyperplane or polytope
    if isa(obj2,'halfspace') || isa(obj2,'conHyperplane') || ...
       isa(obj2,'mptPolytope') || isa(obj2,'ellipsoid')
   
        res = isIntersecting(obj2,obj1,type);
   
    else
        
        % exact or over-approximative algorithm
        if strcmp(type,'exact')
            
            % convert objects to constrained zonotopes
            if isa(obj2,'zonotope') || isa(obj2,'interval')
               obj2 = conZonotope(obj2); 
            end
            
            % conZonotope and conZonotope intersection
            if isa(obj2,'conZonotope')
               res = ~isempty(obj1 & obj2);
            elseif isa(obj2,'zonoBundle')
               res = isIntersecting(obj2,obj1,type); 
            else
               error('No exact algorithm implemented for this set representation!'); 
            end
            
        else       
            if isa(obj2,'inverval')
               res = isIntersecting(obj2,obj1,type); 
            else
               res = isIntersecting(interval(obj1),obj2,type); 
            end
        end
    end
end

%------------- END OF CODE --------------