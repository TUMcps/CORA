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
    
    function han = plot(pgon,dims,varargin)
    % plot the polygon object (always 2D!)
    % input argument dims is here for compliance with other plot functions)
    
        % parse input arguments
        NVpairs = {'EdgeColor',colorblind('b')};
        
        if nargin >= 3
            NVpairs = readPlotOptions(varargin(1:end),'polygon');
            % name-value 'Filled',true|false is deprecated: issue warning
            [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical');
            % faceAlpha = 1 (not the default 0.35 from @polyshape/plot)
            % except otherwise specified by user
            [NVpairs,facealpha] = readNameValuePair(NVpairs,'FaceAlpha','isscalar');
            if ~isempty(filled)
                warning("Name-value pair 'Filled'-true|false is ignored.");
            end
        end
        % readout 'FaceColor' to decide plot/fill call where necessary
        [~,facecolor] = readNameValuePair(NVpairs,'FaceColor');
        
        % plot the polygon object
        if isempty(facecolor) || strcmp(facecolor,'none') || ~isempty(facealpha)
            han = plot(pgon.set,NVpairs{:});            
        else
            han = plot(pgon.set,NVpairs{:},'FaceAlpha',1);
        end
    end
    
    function c = center(pgon)
    % get the center of the polygon
        [x,y] = centroid(pgon.set);
        c = [x;y];
    end
    
    function res = isIntersecting(pgon1,pgon2)
    % check if two polygon object insterst
        temp = intersect(pgon1.set,pgon2.set);
        res = ~isempty(temp.Vertices);
    end
    
    function res = contains(pgon1,pgon2)
    % check if polygon object obj1 contains polygon obj2 or a vector of
    % points (one logical value for each point)
    
        if isnumeric(pgon2)
            
            res = isinterior(pgon1.set,pgon2(1,:),pgon2(2,:))';
            
        elseif isa(pgon2,'polygon')
            
            % compute union
            u = union(pgon1.set,pgon2.set);

            % check if area of obj1 is identical to area of union
            A1 = area(pgon1.set);
            A2 = area(u);

            res = A1 == A2;
        else
            throw(CORAerror('CORA:notSupported',...
                'This set representation is not supported!'));
        end
    end
    
    function pgon = and(pgon1,pgon2)
    % computes the intersection of two polygons
    
        temp = intersect(pgon1.set,pgon2.set);
        V = temp.Vertices;
        pgon = polygon(V(:,1),V(:,2));  
    end
    
    function pgon = or(pgon1,pgon2)
    % computes the union of two polygons
        
        if isempty(pgon1)
            pgon = pgon2; 
        elseif isempty(pgon2)
            pgon = pgon1;
        else
            temp = union(pgon1.set,pgon2.set);
            V = temp.Vertices;
            pgon = polygon(V(:,1),V(:,2));
        end
        
    end
    
    function pgon = plus(pgon,summand)
    % compute the minkowski sum   
        
       % get polygon object
       [pgon,summand] = findClassArg(pgon,summand,'polygon');
       
       % different types of sets
       if isnumeric(summand)
           % translate the polygon
           pgon.set = translate(pgon.set,summand');
       elseif isa(summand,'polygon')
           % compute Minkowski sum
           w = warning();
           warning('off');
           
           T1 = triangulation(pgon);
           T2 = triangulation(summand);
           pgon = [];
           
           for i = 1:length(T1)
               for j = 1:length(T2)
                   V1 = T1{i}.set.Vertices;
                   V2 = T2{j}.set.Vertices;
                   V = [];
                   for k = 1:size(V1,1)
                       V = [V; V2 + V1(k,:)];
                   end
                   pgon = pgon | convHull(polygon(V(:,1),V(:,2)));
               end
           end
           
           warning(w);
       else
           throw(CORAerror('CORA:noops',pgon,summand)); 
       end
    end
    
    function pgon = minus(pgon,minuend)
    % compute the Minkowski difference
    
        % different types of sets
        if isnumeric(minuend)
           pgon.set = translate(pgon.set,-minuend'); 
        elseif isa(minuend,'polygon')
           % compute Mink. diff. by translating the mirrored minuend along
           % the boundary of the polygon
           c = center(minuend);
           minuend = minuend - c;
           pgon = pgon - c;
           temp = minuend.set.Vertices;
           minuend = polygon(-temp(:,1),-temp(:,2));
           V = pgon.set.Vertices;
           V = [V;V(1,:)]';
           diff = [];
           for i = 1:size(V,2)-1
              diff = linComb(minuend + V(:,i),minuend + V(:,i+1)) | diff;
           end
           pgon.set = subtract(pgon.set, diff.set);
        else
            throw(CORAerror('CORA:noops',pgon,minuend));
        end 
    end

    function pgon = minkDiff(pgon,minuend)
    % compute the Minkowski difference
        pgon = minus(pgon,minuend);
    end
    
    function pgon = mtimes(M,pgon)
    % compute linear transformation of a polygon
       
       % check dimension of the matrix
       if ~isnumeric(M) || any(size(M) ~= [2,2])
           throw(CORAerror('CORA:specialError',...
               'Operation "mtimes" is only defined for square matrices of dimension 2!')); 
       end
       
       % multiplication with matrix
       w = warning();
       warning('off');
       
       V = pgon.set.Vertices;
       V = (M*V')'; 

       pgon = polygon(V(:,1),V(:,2));
       
       warning(w);
    end
    
    function pgon = convHull(pgon,varargin)
    % compute convex hull of a polygon

        % only for one other set implemented
        if nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end

        % compute union if two sets are passed
        if nargin > 1
            pgon = pgon | varargin{1};
        end
        
        % compute convex hull
        pgon = convhull(pgon.set);
        
        % construct resulting polygon object
        V = pgon.Vertices;
        pgon = polygon(V(:,1),V(:,2));
    end
    
    function pgon = quadMap(varargin)
    % compute tight enclosure of the quadratic map
    
        if nargin < 2
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin > 3
            throw(CORAerror('CORA:tooManyInputArgs',3));

        elseif nargin == 2
            
            % parse input arguments
            pgon = varargin{1};
            Q = varargin{2};
            
            if ~all(size(Q) == [2,1])
                throw(CORAerror('CORA:wrongValue','second',"be of the size [2,1]"));
            end
            
            % compute triangulation of polygon and compute quadratic map
            % for each triangle using polynomial zonotopes
            list = triangulation(pgon);
            pgon = [];
            
            for i = 1:length(list)
                
               % convert current triangle to polynomial zonotope 
               poly = mptPolytope(list{i}.set.Vertices);
               pZ = polyZonotope(poly);
               
               % compute quadratic map
               temp = quadMap(pZ,Q);
               
               % compute polygon enclousre of the result
               temp = polygon(temp);
               
               % unite with results for other triangles
               pgon = pgon | temp;
               
            end
            
        elseif nargin == 3
            
            % parse input arguments
            pgon1 = varargin{1};
            pgon2 = varargin{2};
            Q = varargin{3};
            
            if ~all(size(Q) == [2,1])
                throw(CORAerror('CORA:wrongValue','second',"be of the size [2,1]"));
            end
            
            % compute triangulation of polygons and compute quadratic map
            % for each triangle combinations using polynomial zonotopes
            list1 = triangulation(pgon1);
            list2 = triangulation(pgon2);
            
            pgon = [];
            
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
                   pgon = pgon | temp;
               
               end
            end
            
        end       
    end
    
    function pgon = linComb(pgon1,pgon2)
    % compute linear combination of two polygons
        
        w = warning();
        warning('off');
    
        list1 = triangulation(pgon1);
        list2 = triangulation(pgon2);
        
        pgon = [];
        
        for i = 1:length(list1)
            for j = 1:length(list2)
                pgon = pgon | convHull(list1{i},list2{j});
            end
        end
        
        warning(w);
    end
        
    function I = interval(pgon)
    % compute an interval enclosure of the polygon
    
        V = pgon.set.Vertices;
        I = interval(min(V,[],1)',max(V,[],1)');
    end

    function P = mptPolytope(pgon)
    % convert a convex polygon to a mptPolytope

        if isConvex(pgon)
            P = mptPolytope(pgon.set.Vertices);
        else
            P = mptPolytope(convHull(pgon));
        end
    end
    
    function p = randPoint(pgon,varargin)
    % draw random points within the polygon
        
        % parse input arguments
        [N,type] = setDefaultValues({1,'standard'},varargin{:});
        
        % different types of random points
        if strcmp(type,'standard')
           
            p = zeros(2,N);
            
            list = triangulation(pgon);
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
            
                p = pgon.set.Vertices';
            
            else
               
                V = pgon.set.Vertices;
                
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
            
        end
    end

    function res = isConvex(pgon)
    % check if the polygon is convex

        res = false;
        try
            % polygon is convex if its convex hull has the same number of
            % vertices as the polygon itself already has
            if length(convhull(pgon.set.Vertices))-1 == size(pgon.set.Vertices,1)
                res = true;
            end
        catch
            warning("Built-in method 'convHull' failed. Result returned as 'false'.");
        end
    end

    function res = dim(pgon)
    % return dimension of the polygon (always two-dimensional)

        res = 2;
    end
    
    function list = triangulation(pgon)
    % compute an triangulation of the polygon
    
        % compute triangulation
        T = triangulation(pgon.set);
        
        % convert the triangles to polygons
        list = cell(size(T.ConnectivityList,1),1);
        
        for i = 1:size(T.ConnectivityList,1)
           V = T.Points(T.ConnectivityList(i,:),:);
           list{i} = polygon(V(:,1),V(:,2));
        end
    end
    
    function pgon = simplify(pgon,varargin)
    % enclose the polygon by a simpler polygon with less vertices   
        
        if nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end
        % default values
        tol = setDefaultValues({0.01},varargin{:});
        % check input arguments
        inputArgsCheck({{pgon,'att','polygon'},...
                        {tol,'att','numeric',{'scalar','nonnegative'}}});
    
        % simplify polygon boundary using the Douglar-Peucker algorithm 
        V = pgon.set.Vertices;
        V_ = douglasPeucker(V',tol);
        pgon = polygon(V_(1,:),V_(2,:));
        
        % enlarge the simplified polygon so that it is guaranteed that it
        % encloses the original polygon
        temp = polygon([-tol -tol tol tol],[-tol tol tol -tol]);
        pgon = pgon + temp;
    end

    function list = splitIntoConvexSets(pgon)
    % split a polygon into a set of convex regions

        % compute triangulation
        list = triangulation(pgon);

        % unite neighbouring polygons to obtain larger convex shapes
        while true
    
            finished = true;

            % loop over all polygons in the current list
            for i = 1:length(list)

                stop = false;

                % try to combine it with other polygons in the list
                for j = i+1:length(list)

                    tmp = list{i} | list{j};

                    if isConvex(tmp)
                        list{i} = tmp;
                        list{j} = [];
                        stop = true; finished = false;
                        break;
                    end
                end

                if stop
                    break;
                end
            end

            % remove empty list entries
            list = list(~cellfun('isempty',list));

            if finished
                break;
            end
        end
    end
end

methods (Static = true)
    
    function pgon = generateRandom() 
    % generate random polygon
        
        points = 10*rand(1) * (-1 + 2*rand(2,100)) + 10*(-1 + 2*rand(2,1));
        pgon = polygon.enclosePoints(points);
    end
    
    function pgon = enclosePoints(points) 
    % enclose point cloud with polygon
    
        ind = boundary(points',0.5);
        pgon = polygon(points(1,ind),points(2,ind));
    end
end

end

%------------- END OF CODE --------------