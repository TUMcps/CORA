function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if interval obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - interval object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    int1 = interval([0;0],[2;2]);
%    int2 = interval([1;1],[3;3]);
%    int3 = interval([-3;-3],[-1;1]);
%
%    isIntersecting(int1,int2)
%    isIntersecting(int1,int3)
%
%    figure
%    hold on
%    plot(int1,[1,2],'b');
%    plot(int2,[1,2],'g');
%
%    figure
%    hold on
%    plot(int1,[1,2],'b');
%    plot(int3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isIntersecting

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      22-July-2016
% Last update:  14-Sep-2019
%               21-Nov-2019 (NK, added intersection with other sets)
%               01-July-2021 (MW, bug fix in 1D intersection)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && ~isempty(varargin{1}) 
        type = varargin{1};
    end
    
    % get interval object
    if ~isa(obj1,'interval')
       temp = obj1;
       obj1 = obj2;
       obj2 = temp;
    end
    
    % interval and interval intersection
    if isa(obj2,'interval')
        
        res = true;
        
        % get object properties
        sup1 = obj1.sup;
        inf1 = obj1.inf;
        sup2 = obj2.sup;
        inf2 = obj2.inf;
        
        % loop over all dimensions
        for i = 1:length(obj1)
           if ~isIntersecting1D(inf1(i),sup1(i),inf2(i),sup2(i))
              res = false;
              return
           end
        end
        
    elseif isa(obj2,'halfspace') || isa(obj2,'conHyperplane') || ...
           isa(obj2,'mptPolytope')
        
        res = isIntersecting(obj2,obj1,type);
        
    else
        
        % exact or over-approximative algorithm
        if strcmp(type,'exact')
           
            res = isIntersecting(obj2,obj1,type);
            
        else
            
            res = isIntersecting(mptPolytope(obj1),obj2,type);
            
        end
    end
    
end


% Auxiliary Functions -----------------------------------------------------

function res = isIntersecting1D(inf1,sup1,inf2,sup2)
% check if two one-dimensional intervals intersect

    res = false;

    if inf1 <= inf2
        if inf2 <= sup1
            res = true;
        end
        
    else % inf2 < inf1
        if inf1 <= sup2
            res = true;
        end
        
    end
    
end

%------------- END OF CODE --------------