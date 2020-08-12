classdef polygon
% polygon class 
%
% Syntax:  
%    obj = polygon(x,y)
%
% Inputs:
%    x - vector with x coordinates of the polygon vertices
%    y - vector with y coordinates of the polygon vertices
%
% Outputs:
%    obj - polygon object
%
% Example:
%    x = gallery('uniformdata',30,1,1);
%    y = gallery('uniformdata',30,1,10);
%    ind = boundary(x,y);
%
%    pgon = polygon(x(ind),y(ind));
%
%    plot(pgon,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Niklas Kochdumper
% Written:      13-March-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    set = []
end
    
methods
    
    % class constructor
    function obj = polygon(x,y)
        obj.set = polyshape(x,y);
    end
    
    function han = plot(obj,dims,varargin)
    % plot the polygon object (always 2D!)
    % input argument dims is here for compliance with other plot functions)
    
        % parse input arguments
        filled = 0;
        linespec = 'b';
        NVpairs = {};
        
        if nargin >= 3
            [linespec,NVpairs] = readPlotOptions(varargin(1:end));
            [NVpairs,filled] = readNameValuePair(NVpairs,'Filled', ...
                                                 'islogical',filled);
        end
        
        % plot the polygon object
        if filled
            han = plot(obj.set,'FaceColor',linespec,'FaceAlpha',...
                       1,NVpairs{:});
        else
            han = plot(obj.set,'FaceColor','none','EdgeColor', ...
                       linespec,NVpairs{:});
        end
    end
    
    function c = center(obj)
    % get the center of the polygon
        [x,y] = centroid(obj.set);
        c = [x;y];
    end
    
    function res = isIntersecting(obj1,obj2)
    % check if two polygon object insterst
        temp = intersect(obj1.set,obj2.set);
        res = ~isempty(temp.Vertices);
    end
    
    function res = in(obj1,obj2)
    % check if polygon obj2 is inside polygon obj2
    
        if isnumeric(obj2)
            
            res = isinterior(obj1.set,obj2(1),obj2(2));
            
        elseif isa(obj2,'polygon')
            
            % compute union
            u = union(obj1.set,obj2.set);

            % check if area of obj1 is identical to area of union
            A1 = area(obj1.set);
            A2 = area(u);

            res = A1 == A2;
        else
            error('This set representation is not supported!');
        end
    end
    
    function obj = or(obj1,obj2)
    % computes the union of two polygons
        
        if isempty(obj1)
            obj = obj2; 
        elseif isempty(obj2)
            obj = obj1;
        else
            temp = union(obj1.set,obj2.set);
            V = temp.Vertices;
            obj = polygon(V(:,1),V(:,2));
        end
        
    end
    
    function obj = plus(obj,summand)
    % compute the minkowski sum with a point    
        
       % get polygon object
       if ~isa(obj,'polygon')
          temp = obj;
          obj = summand;
          summand = temp;
       end
       
       % check if second argument is a vector
       if ~isnumeric(summand)
          error('Operation "plus" is not yet implemented for this set representation!'); 
       end
       
       % translate the polygon
       obj.set = translate(obj.set,summand');
    end
    
    function obj = mtimes(mat,obj)
    % compute linear transformation of a polgon
       
       % check dimension of the matrix
       if ~isnumeric(mat) || size(mat,1) ~= 2 || size(mat,2) ~= 2
          error('Operation "mtimes" is only defined for square matrices of dimension 2!'); 
       end
       
       % multiplication with matrix
       w = warning();
       warning('off');
       
       V = obj.set.Vertices;
       V = (mat*V')'; 

       obj = polygon(V(:,1),V(:,2));
       
       warning(w);
    end
    
    function int = interval(obj)
    % compute an interval enclosure of the polygon
    
        V = obj.set.Vertices;
        int = interval(min(V,[],1)',max(V,[],1)');
    end
end
end

%------------- END OF CODE --------------