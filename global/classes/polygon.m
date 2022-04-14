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
    
    function obj = and(obj1,obj2)
    % computes the intersection of two polygons
    
        temp = intersect(obj1.set,obj2.set);
        V = temp.Vertices;
        obj = polygon(V(:,1),V(:,2));  
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
       
       % different types of sets
       if isnumeric(summand)
           % translate the polygon
           obj.set = translate(obj.set,summand');
       elseif isa(summand,'polygon')
           % compute Minkowski sum
           w = warning();
           warning('off');
           
           T1 = triangulation(obj);
           T2 = triangulation(summand);
           obj = [];
           
           for i = 1:length(T1)
               for j = 1:length(T2)
                   V1 = T1{i}.set.Vertices;
                   V2 = T2{j}.set.Vertices;
                   V = [];
                   for k = 1:size(V1,1)
                       V = [V; V2 + V1(k,:)];
                   end
                   obj = obj | convHull(polygon(V(:,1),V(:,2)));
               end
           end
           
           warning(w);
       else
          error('Operation "plus" is not yet implemented for this set representation!'); 
       end
    end
    
    function obj = mtimes(mat,obj)
    % compute linear transformation of a polygon
       
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
    
    function obj = convHull(obj1,varargin)
    % compute convex hull of a polygon
    
        % compute union if two sets are passed
        if nargin > 1
            temp = obj1 | varargin{1};
        else
            temp = obj1; 
        end
        
        % compute convex hull
        temp = convhull(temp.set);
        
        % construct resulting polygon object
        V = temp.Vertices;
        obj = polygon(V(:,1),V(:,2));
    end
    
    function obj = quadMap(varargin)
    % compute tight enclosure of the quadratic map
    
        if nargin == 2
            
            % parse input arguments
            pgon = varargin{1};
            Q = varargin{2};
            
            if ~all(size(Q) == [2,1])
               [msg,id] = errWrongInput('Q');
               error(msg,id);
            end
            
            % compute triangulation of polygon and compute quadratic map
            % for each triangle using polynomial zonotopes
            list = triangulation(pgon);
            obj = [];
            
            for i = 1:length(list)
                
               % convert current triangle to polynomial zonotope 
               poly = mptPolytope(list{i}.set.Vertices);
               pZ = polyZonotope(poly);
               
               % compute quadratic map
               temp = quadMap(pZ,Q);
               
               % compute polygon enclousre of the result
               temp = polygon(temp);
               
               % unite with results for other triangles
               obj = obj | temp;
               
            end
            
        elseif nargin == 3
            
            % parse input arguments
            pgon1 = varargin{1};
            pgon2 = varargin{2};
            Q = varargin{3};
            
            if ~all(size(Q) == [2,1])
               [msg,id] = errWrongInput('Q');
               error(msg,id);
            end
            
            % compute triangulation of polygons and compute quadratic map
            % for each triangle combinations using polynomial zonotopes
            list1 = triangulation(pgon1);
            list2 = triangulation(pgon2);
            
            obj = [];
            
            for i = 1:length(list1)
                
               % convert current triangle to polynomial zonotope 
               poly = mptPolytope(list1{i}.set.Vertices);
               pZ1 = polyZonotope(poly);
               
               for j = 1:length(list2)
               
                   % convert current triangle to polynomial zonotope 
                   poly = mptPolytope(list2{i}.set.Vertices);
                   pZ2 = polyZonotope(poly);
                   
                   % compute quadratic map
                   temp = quadMap(pZ1,pZ2,Q);

                   % compute polygon enclosure of the result
                   temp = polygon(temp);

                   % unite with results for other triangles
                   obj = obj | temp;
               
               end
            end
            
        else
            error('Wrong number of input arguments')
        end       
    end
    
    function obj = linComb(obj1,obj2)
    % compute linear combination of two polygons
        
        w = warning();
        warning('off');
    
        list1 = triangulation(obj1);
        list2 = triangulation(obj2);
        
        obj = [];
        
        for i = 1:length(list1)
            for j = 1:length(list2)
                obj = obj | convHull(list1{i},list2{j});
            end
        end
        
        warning(w);
    end
        
    function int = interval(obj)
    % compute an interval enclosure of the polygon
    
        V = obj.set.Vertices;
        int = interval(min(V,[],1)',max(V,[],1)');
    end
    
    function p = randPoint(obj,varargin)
    % draw random points within the polygon
        
        % parse input arguments
        N = 1;
        type = 'normal';
        if nargin > 1 && ~isempty(varargin{1})
           N = varargin{1}; 
        end
        if nargin > 2 && ~isempty(varargin{2})
           type = varargin{2}; 
        end
        
        % different types of random points
        if strcmp(type,'normal')
           
            p = zeros(2,N);
            
            list = triangulation(obj);
            cnt = 1;
            cnt_ = 1;
            
            while cnt <= N
                
                % get random point by interpolation between the vertices
                V = list{cnt_}.set.Vertices;
                d = rand(3,1); d = d./sum(d);
                p(:,cnt) = sum(V' * diag(d),2);
      
                % update counter
                cnt = cnt + 1;
                cnt_ = cnt_ + 1;
                if cnt_ > length(list)
                   cnt_ = 1; 
                end
            end
            
        elseif strcmp(type,'extreme')
           
            % return all extreme point
            if ischar(N) && strcmp(N,'all')
            
                p = obj.set.Vertices';
            
            else
               
                V = obj.set.Vertices;
                
                if N <= size(V,1)
                    ind = randperm(size(V,1));
                    ind = ind(1:N);
                    p = V(ind,:)';
                else
                    p = [V',zeros(2,size(V,1)-N)];
                    for i = size(V,1)+1:N
                        ind = randi(size(V,1)-1,[1,1]);
                        d = rand(2,1); d = d/sum(d);
                        v1 = V(ind,:)'; v2 = V(ind+1,:)';
                        p(:,i) = d(1)*v1 + d(2)*v2;
                    end
                end
            end
            
        else
            [msg,id] = errWrongInput('type');
            error(id,msg);
        end
    end
    
    function list = triangulation(obj)
    % compute an triangulation of the polygon
    
        % compute triangulation
        T = triangulation(obj.set);
        
        % convert the triangles to polygons
        list = cell(size(T.ConnectivityList,1),1);
        
        for i = 1:size(T.ConnectivityList,1)
           V = T.Points(T.ConnectivityList(i,:),:);
           list{i} = polygon(V(:,1),V(:,2));
        end
    end
    
    function obj = simplify(obj,varargin)
    % enclose the polygon by a simpler polygon with less vertices   
        
        % parse input arguments
        tol = 0.01;
        if nargin >= 2 && ~isempty(varargin{1})
            tol = varargin{1};
        end
    
        % simplify polygon boundary using the Douglar-Peucker algorithm 
        V = obj.set.Vertices;
        V_ = douglasPeucker(V',tol);
        obj = polygon(V_(1,:),V_(2,:));
        
        % enlare the simplified polygon so that it is guaranteed that it
        % encloses the original polygon
        temp = polygon([-tol -tol tol tol],[-tol tol tol -tol]);
        obj = obj + temp;
    end
end

methods (Static = true)
    
    function obj = generateRandom() 
    % generate random polygon
        
        points = 10*rand(1) * (-1 + 2*rand(2,100)) + 10*(-1 + 2*rand(2,1));
        obj = polygon.enclosePoints(points);
    end
    
    function obj = enclosePoints(points) 
    % enclose point cloud with zonotope
    
        ind = boundary(points',0.5);
        obj = polygon(points(1,ind),points(2,ind));
    end
end
end

%------------- END OF CODE --------------