function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if polyZonotope obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - polyZonotope object
%    obj2 - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    pZ = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    hs1 = halfspace([-1 1],-2);
%    hs2 = halfspace([-1 1],-4);
%
%    isIntersecting(pZ,hs1,'approx')
%    isIntersecting(pZ,hs2,'approx')
%
%    figure
%    hold on
%    xlim([-4,6]);
%    ylim([-5,5]);
%    plot(hs1,[1,2],'b');
%    plot(pZ,[1,2],'g','Filled',true,'EdgeColor','none');
%
%    figure
%    hold on
%    xlim([-4,6]);
%    ylim([-5,5]);
%    plot(hs2,[1,2],'b');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse optional input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    if strcmp(type,'exact')
       error('No exact algorithm implemented for polyZonotope object!'); 
    end

    % get polyZonotope object
    if ~isa(obj1,'polyZonotope')
        temp = obj1;
        obj1 = obj2;
        obj2 = temp;
    end
    
    % call function for other set representations
    if isa(obj2,'polyZonotope')
        
        res = isIntersecting(zonotope(obj1),zonotope(obj2),'approx');
        
    elseif isa(obj2,'halfspace') || isa(obj2,'conHyperplane') || ...
           isa(obj2,'mptPolytope') || isa(obj2,'zonotope') || ...
           isa(obj2,'interval') || isa(obj2,'zonoBundle') || ...
           isa(obj2,'conZonotope')
   
        res = isIntersecting(obj2,obj1,type);
   
    else
        
        error('Operation not implemented!');
        
    end     
end

%------------- END OF CODE --------------