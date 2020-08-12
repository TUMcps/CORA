function res = in(Z,obj,varargin)
% in - determines if obj is contained in Z
%
% Syntax:  
%    res = in(Z,obj)
%    res = in(Z,obj,type)
%
% Inputs:
%    Z - zonotope object
%    obj - contSet object or point
%    type - type of containment check ('exact' or 'approx')
%
% Outputs:
%    res - boolean whether obj is contained in Z, or not
%
% Example: 
%    zono1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    zono2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    zono3 = zono2 + [3;0];
% 
%    in(zono1,zono2)
%    in(zono1,zono3)
% 
%    figure
%    hold on
%    plot(zono1,[1,2],'b');
%    plot(zono2,[1,2],'g');
%    
%    figure
%    hold on
%    plot(zono1,[1,2],'b');
%    plot(zono3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/in, conZonotope/in
%
% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      07-May-2007 
% Last update:  06-April-2017
%               14-Sep-2019
%               19-Nov-2019 (NK, changed to header format)
% Last revision:---

%------------- BEGIN CODE --------------

    res = true;

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && strcmp(varargin{1},'approx')
        type = 'approx';
    end
        
    % point in zonotope containment
    if isnumeric(obj)
        
        Z = conZonotope(Z);
        
        for i = 1:size(obj,2)
            res = in(conZonotope(Z),obj(:,i));
            if res ~= 1
               return; 
            end
        end
        
    % capsule/ellipsoid in zonotope containment
    elseif isa(obj,'capsule') || isa(obj,'ellipsoid')
        
        poly = mptPolytope(Z);
        res = in(poly,obj); 

    else
        
        % use the fast but over-approximative or the exact but possibly
        % slow containment check
        if strcmp(type,'exact')

            if isa(obj,'taylm') || isa(obj,'polyZonotope')
                error('Exact containment check not possible for this set representation!'); 
            elseif isa(obj,'interval')
                res = in(Z,vertices(obj));
            else
                poly = mptPolytope(Z);
                res = in(poly,obj); 
            end
            
        else
            
            if isa(obj,'taylm') || isa(obj,'polyZonotope')
                poly = mptPolytope(Z);
                res = in(poly,obj); 
            else
                cZ = conZonotope(Z);
                res = in(cZ,obj);
            end
        end
    end
end

%------------- END OF CODE --------------