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
%    c = [0;0];
%    G = [8 4; 0 8];
%    expMat = [1 0; 0 1; 0 0];
%    A = [1 1 -0.25];
%    b = 0.75;
%    expMat_ = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%    pZ1 = polyZonotope([0;0],0.5*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%
%    res1 = isIntersecting(cPZ,pZ1,'approx')
%    res2 = isIntersecting(cPZ,pZ2,'approx')
%
%    figure; hold on;
%    plot(cPZ,[1,2],'b','Filled',true,'EdgeColor','none','Splits',15);
%    plot(pZ1,[1,2],'g');
%    plot(pZ2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/isIntersecting, isempty

% Author:       Niklas Kochdumper
% Written:      04-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse optional input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    if strcmp(type,'exact')
       error(errNoExactAlg(obj2,obj1));
    end

    % get polyZonotope object
    if ~isa(obj1,'conPolyZono')
        temp = obj1; obj1 = obj2; obj2 = temp;
    end
    
    % call function for other set representations
    if isa(obj2,'conPolyZono') || isa(obj2,'polyZonotope') || ...
       isa(obj2,'capsule')
        
        % fast check based on zonotope enclosure
        res = isIntersecting(zonotope(obj1),zonotope(obj2));
        
        if ~res
            return; 
        end
        
        % convert second set to constrained polynomial zonotope
        obj2 = conPolyZono(obj2);
        
        % compute intersection of the two sets
        int = obj1 & obj2;
        
        % check if the intersection is empty
        res = ~isempty(int);
        
        
    elseif isa(obj2,'halfspace') || isa(obj2,'conHyperplane') || ...
           isa(obj2,'mptPolytope') || isa(obj2,'zonotope') || ...
           isa(obj2,'interval') || isa(obj2,'zonoBundle') || ...
           isa(obj2,'conZonotope') || isa(obj2,'ellipsoid')
   
        res = isIntersecting(obj2,obj1,type);
   
    else
        
        % throw error for given arguments
        error(noops(obj1,obj2));
    end     
end

%------------- END OF CODE --------------