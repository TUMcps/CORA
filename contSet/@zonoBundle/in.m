function res = in(obj1,obj2,varargin)
% in - determines if obj2 is contained in obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%    res = in(obj1,obj2,type)
%
% Inputs:
%    obj1 - zonoBundle object
%    obj2 - contSet object or point
%    type - type of containment check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if obj2 is contained in obj1, or not
%
% Example: 
%    int1 = interval([0;-1],[2;2]);
%    int2 = int1 + [2;0];
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    in(zB,int1)
%    in(zB,int2)
%
%    figure
%    hold on
%    plot(zB,[1,2],'b');
%    plot(int1,[1,2],'g');
%    
%    figure
%    hold on
%    plot(zB,[1,2],'b');
%    plot(int2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/in, conZonotope/in
%
% Author:       Niklas Kochdumper
% Written:      19-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && strcmp(varargin{1},'approx')
        type = 'approx';
    end
        
    % point or point cloud in zonotope bundle containment
    if isnumeric(obj2)
        
        for i = 1:size(obj2,2)
            res = in(conZonotope(obj1),obj2(:,i));
            if res ~= 1
                return;
            end
        end
        
    % capsule/ellipsoid in zonotope bundle containment
    elseif isa(obj2,'capsule') || isa(obj2,'ellipsoid')
        
        poly = mptPolytope(obj1);
        res = in(poly,obj2); 

    else
        
        % use the fast but over-approximative or the exact but possibly
        % slow containment check
        if strcmp(type,'exact')

            if isa(obj2,'taylm') || isa(obj2,'polyZonotope')
                error('Exact containment check not possible for this set representation!'); 
            elseif isa(obj2,'interval')
                res = in(obj1,vertices(obj2));
            else
                poly = mptPolytope(obj1);
                res = in(poly,obj2); 
            end
            
        else
            
            if isa(obj2,'taylm') || isa(obj2,'polyZonotope')
                poly = mptPolytope(obj1);
                res = in(poly,obj2); 
            else
                cZ = conZonotope(obj1);
                res = in(cZ,obj2,type);
            end
        end
    end
end

%------------- END OF CODE --------------