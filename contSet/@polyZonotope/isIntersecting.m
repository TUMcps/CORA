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
%    pZ1 = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    pZ2 = polyZonotope([1;4],[1 1 1 1;1 0 -1 1],[],[1 0 1 2;0 1 1 0]);
%    hs1 = halfspace([-1 1],-2);
%    hs2 = halfspace([-1 1],-4);
%
%    isIntersecting(pZ1,pZ2,'approx')
%    isIntersecting(pZ1,hs1,'approx')
%    isIntersecting(pZ1,hs2,'approx')
%
%    figure; hold on;
%    plot(pZ1,[1,2],'b');
%    plot(pZ2,[1,2],'r');
%
%    figure; hold on
%    xlim([-4,6]); ylim([-5,5]);
%    plot(hs1,[1,2],'b');
%    plot(pZ1,[1,2],'g','Filled',true,'EdgeColor','none');
%
%    figure; hold on
%    xlim([-4,6]); ylim([-5,5]);
%    plot(hs2,[1,2],'b');
%    plot(pZ1,[1,2],'r','Filled',true,'EdgeColor','none');
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
       error(errNoExactAlg(obj2,obj1));
    end

    % get polyZonotope object
    if ~isa(obj1,'polyZonotope')
        temp = obj1;
        obj1 = obj2;
        obj2 = temp;
    end
    
    % call function for other set representations
    if isa(obj2,'polyZonotope')
        
        % fast check based on zonotope enclosure
        res = isIntersecting(zonotope(obj1),zonotope(obj2));
        
        if ~res
            return; 
        end
        
        % construct polynomial constraint for the intersection
        c = obj1.c - obj2.c;
        G = [obj1.G -obj2.G];
        Grest = [obj1.Grest,-obj2.Grest];
        expMat = blkdiag(obj1.expMat,obj2.expMat);
        
        % contract the factor domain \alpha_k in [-1,1] based on the
        % polynomial constraint
        n = size(expMat,1) + size(Grest,2);
        dom = interval(-ones(n,1),ones(n,1));
        
        int = contractPoly(c,G,Grest,expMat,dom,'linearize',3,7);
        
        % check if polynomial zonotopes are intersecting
        if isempty(int)
           res = false;
        else
           res = true;
        end
        
    elseif isa(obj2,'halfspace') || isa(obj2,'conHyperplane') || ...
           isa(obj2,'mptPolytope') || isa(obj2,'zonotope') || ...
           isa(obj2,'interval') || isa(obj2,'zonoBundle') || ...
           isa(obj2,'conZonotope') || isa(obj2,'ellipsoid') || ...
           isa(obj2,'conPolyZono')
   
        res = isIntersecting(obj2,obj1,type);
   
    else
        
        % throw error for given arguments
        error(noops(obj1,obj2));
    end     
end

%------------- END OF CODE --------------
