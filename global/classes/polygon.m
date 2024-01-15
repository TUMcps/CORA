classdef polygon
% polygon class
%
% Syntax:
%    obj = polygon(x,y)
%    obj = polygon(V)
%    obj = polygon(set)
%
% Inputs:
%    x - vector with x coordinates of the polygon vertices
%    y - vector with y coordinates of the polygon vertices
%    V - vertices (2-dimensional)
%    set - polyshape object
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
% See also: interval, polytope, polyshape

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   09-May-2023 (TL, constructor, plotPolygon)
%                28-June-2023 (TL, minkDiff, minus, uplus, uminus)
%                27-May-2023 (MW, isequal, removeCollinearPoints)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    set = []
end

methods

    % class constructor
    function obj = polygon(varargin)
        if nargin == 0
            obj.set = polyshape();
        elseif nargin == 1
            set = varargin{1};
            inputArgsCheck({{set, 'att', {'polygon', 'polyshape', 'numeric'}}})
            if isa(set, 'polygon')
                obj = set;
            elseif isa(set, 'polyshape')
                obj.set = set;
            else
                V = set;
                if (size(V, 1) ~= 2)
                    throw(CORAerror("CORA:wrongValue", ...
                        'Given vertices should be two dimensional.'))
                end
                obj.set = polyshape(V(1, :), V(2, :));
            end
        elseif nargin == 2
            x = varargin{1};
            y = varargin{2};
            inputArgsCheck({{x, 'att', 'numeric'}, {y, 'att', 'numeric'}})
            if ~isvector(x) || ~all(size(x) == size(y))
                throw(CORAerror("CORA:wrongInputInConstructor", ...
                    'Given vertices x,y need to be vectors of the same length.'))
            end
            obj.set = polyshape(x, y);
        else
            throw(CORAerror('CORA:tooManyInputArgs'));
        end
    end

    function han = plot(pgon, dims, varargin)
        % plot the polygon object (always 2D!)
        % input argument dims is here for compliance with other plot functions)

        if nargin < 2
            dims = [1, 2];
        end
        n = length(dims);

        if pgon.set.NumRegions > 2
            % plot each region separately
            regions = pgon.set.regions;
            pgons = cell(length(regions),1);
            for i=1:length(regions)
                pgons{i} = regions(i);
            end

            han = plotMultipleSetsAsOne(pgons, dims, varargin);

        elseif pgon.set.NumRegions == 1
            % init region
            regions = pgon.set.regions;
            region = regions(1);
        
            % save holes
            holes = region.holes;
            
            % remove holes
            region = region.rmholes();

            % get polygon of region
            V = [region.Vertices;region.Vertices(1,:)]';

            % cut out holes; inserted at smallest distance to boundary
            for i=1:length(holes)
                % extract
                holesV = [holes(i).Vertices;holes(i).Vertices(1,:)]';
                holesV3D = reshape(holesV,n,1,[]);

                [~,idxV] = min(min(vecnorm(V-holesV3D,2,1),[],3),[],2);
                [~,idxH] = min(min(vecnorm(V-holesV3D,2,1),[],2),[],3);

                V = [V(:,1:idxV), ...
                    holesV(:,idxH:end),holesV(:,1:idxH), ...
                    V(:,idxV:end)];
            end 

            % plot
            han = plotPolygon(V, varargin{:});

        else
            % plot empty
            V = zeros(2, 0);
            han = plotPolygon(V, varargin{:});
        end

        if nargout == 0
            clear han
        end
    end

    function c = center(pgon)
        % get the center of the polygon
        [x, y] = centroid(pgon.set);
        c = [x; y];
    end

    function res = isIntersecting(pgon1, pgon2)
        % check if two polygon objects intersect
        temp = intersect(pgon1.set, pgon2.set);
        res = ~isempty(temp.Vertices);
    end

    function res = contains(pgon1, pgon2)
        % check if polygon object obj1 contains polygon obj2 or a vector of
        % points (one logical value for each point)

        if isnumeric(pgon2)

            res = isinterior(pgon1.set, pgon2(1, :), pgon2(2, :))';

        elseif isa(pgon2, 'polygon')

            % compute union
            u = union(pgon1.set, pgon2.set);

            % check if area of obj1 is identical to area of union
            A1 = area(pgon1.set);
            A2 = area(u);

            res = A1 == A2;
        else
            throw(CORAerror('CORA:notSupported', ...
                'This set representation is not supported!'));
        end
    end

    function pgon = and(pgon1, pgon2)
        % computes the intersection of two polygons

        temp = intersect(pgon1.set, pgon2.set);
        V = temp.Vertices;
        pgon = polygon(V(:, 1), V(:, 2));
    end

    function pgon = or(pgon1, pgon2)
        % computes the union of two polygons

        if representsa_(pgon1,'emptySet',eps)
            pgon = pgon2;
        elseif representsa_(pgon2,'emptySet',eps)
            pgon = pgon1;
        else
            temp = union(pgon1.set, pgon2.set);
            V = temp.Vertices;
            pgon = polygon(V(:, 1), V(:, 2));
        end

    end

    function pgon = plus(pgon, summand)
        % compute the minkowski sum

        % get polygon object
        [pgon, summand] = findClassArg(pgon, summand, 'polygon');

        % different types of sets
        if isnumeric(summand)
            if isscalar(summand)
                % change to vector
                summand = [summand;summand];
            end

            % translate the polygon
            pgon.set = translate(pgon.set, summand');

        elseif isa(summand, 'polygon')
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
                    for k = 1:size(V1, 1)
                        V = [V; V2 + V1(k, :)];
                    end
                    pgon = pgon | convHull(polygon(V(:, 1), V(:, 2)));
                end
            end

            warning(w);
        else
            throw(CORAerror('CORA:noops', pgon, summand));
        end
    end

    function pgon = minus(pgon, subtrahend)
        % translation of a polygon P by a vector v

        % different types of sets
        if isa(pgon,'polygon') && isnumeric(subtrahend)
            % P - v
            pgon = pgon + -subtrahend;

        elseif isnumeric(pgon) && isa(subtrahend,'polygon')
           % v - P
           pgon = pgon + -subtrahend;

        elseif isa(pgon, 'polygon') && isa(subtrahend, 'polygon')
            % P1 - P2
            % throw error
            classname = class(S);
            throw(CORAerror('CORA:notSupported',...
                sprintf(['The function ''minus'' is not implemented for the class %s except for vectors as a subtrahend.\n', ...
                'If you require to compute the Minkowski difference, use ''minkDiff'' instead.'],classname)));
           
        else
            throw(CORAerror('CORA:noops', pgon, subtrahend));
        end
    end

    function pgon = minkDiff(pgon, subtrahend)
        % compute the Minkowski difference

        % parse input
        if nargin < 2
            throw(CORAerror("CORA:notEnoughInputArgs",2))
        elseif nargin > 2
            throw(CORAerror("CORA:tooManyInputArgs",2))
        end

        inputArgsCheck({ ...
            {pgon,'att','polygon'}; ...
            {subtrahend,'att',{'polygon','numeric'}}
        })

        % numeric case
        if isnumeric(subtrahend)
            pgon = pgon - subtrahend;
            return
        end
        
        % compute Mink. diff. by translating the mirrored subtrahend along
        % the boundary of the polygon
        
        % shift by center
        c = center(subtrahend);
        subtrahend = subtrahend - c;
        pgon = pgon - c;
        
        % mirror subtrahend
        subtrahend = -subtrahend;

        % get boundary of pgon
        V = pgon.set.Vertices';
        
        % translate subtrahend along boundary
        diff = [];
        boundStart = 1;
        vs = size(V, 2);
        for i = 1:vs
            % boundaries are splitted by nan values
            
            % take first point on boundary
            Vi = V(:, i);
            if all(isnan(Vi))
                % region boundary
                continue
            end
        
            % get subsequent point on boundary
            if i+1 <= vs
                % take next value
                Vi1 = V(:, i+1);
            else
                % take start point of current boundary
                Vi1 = V(:, boundStart);
            end
        
            % check if end of current boundary is reached
            if all(isnan(Vi1))
                % take start point of current boundary
                Vi1 = V(:, boundStart);
                % shift to start of next boundary
                boundStart = i + 2;
            end
        
            % linear combination between the two points
            diff = linComb(subtrahend+Vi, subtrahend+Vi1) | diff;
        end
        pgon.set = subtract(pgon.set, diff.set);
    end

    function pgon = mtimes(M, pgon)
        % compute linear transformation of a polygon

        % check dimension of the matrix
        if ~isnumeric(M) && ~(all(size(M) == [2, 2]) || isscalar(M))
            throw(CORAerror('CORA:notSupported', ...
                'Operation "mtimes" is only defined for square matrices of dimension 2 or scalars!'));
        end

        % multiplication with matrix
        w = warning();
        warning('off');

        V = pgon.set.Vertices;
        V = (M * V')';

        pgon = polygon(V(:, 1), V(:, 2));

        warning(w);
    end

    function pgon = times(pgon,M)

        % find class argument
        [pgon,M] = findClassArg(pgon,M,'polygon');

        % check dimension of the matrix
        if ~isnumeric(M) && ~(all(size(M) == [2, 1]) || isscalar(M))
            throw(CORAerror('CORA:notSupported', ...
                'Operation "times" is only defined for vectors matrices of dimension 2 or scalars!'));
        end

        % use mtimes
        pgon = diag(M) * pgon;
    end

    function pgon = uplus(pgon)
        % unary plus
        % pgon = pgon; % identity
    end

    function pgon = uminus(pgon)
        % unary minus
        pgon = -1 * pgon;
    end

    function pgon = convHull(pgon, varargin)
        % compute convex hull of a polygon

        % only for one other set implemented
        if nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs', 2));
        end

        % compute union if two sets are passed
        if nargin > 1
            pgon = pgon | varargin{1};
        end

        % compute convex hull
        pgon = convhull(pgon.set);

        % construct resulting polygon object
        V = pgon.Vertices;
        pgon = polygon(V(:, 1), V(:, 2));
    end

    function pgon = quadMap(varargin)
        % compute tight enclosure of the quadratic map

        if nargin < 2
            throw(CORAerror('CORA:notEnoughInputArgs', 2));
        elseif nargin > 3
            throw(CORAerror('CORA:tooManyInputArgs', 3));

        elseif nargin == 2

            % parse input arguments
            pgon = varargin{1};
            Q = varargin{2};

            if ~all(size(Q) == [2, 1])
                throw(CORAerror('CORA:wrongValue', 'second', "be of the size [2,1]"));
            end

            % compute triangulation of polygon and compute quadratic map
            % for each triangle using polynomial zonotopes
            list = triangulation(pgon);
            pgon = [];

            for i = 1:length(list)
                % convert current triangle to polynomial zonotope 
                poly = polytope(list{i}.set.Vertices');
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

            if ~all(size(Q) == [2, 1])
                throw(CORAerror('CORA:wrongValue', 'second', "be of the size [2,1]"));
            end

            % compute triangulation of polygons and compute quadratic map
            % for each triangle combinations using polynomial zonotopes
            list1 = triangulation(pgon1);
            list2 = triangulation(pgon2);

            pgon = [];

            for i = 1:length(list1)
                
                % convert current triangle to polynomial zonotope 
                poly = polytope(list1{i}.set.Vertices);
                pZ1 = polyZonotope(poly);
                
                for j = 1:length(list2)
               
                    % convert current triangle to polynomial zonotope 
                    poly = polytope(list2{i}.set.Vertices);
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

    function pgon = linComb(pgon1, pgon2)
        % compute linear combination of two polygons

        w = warning();
        warning('off');

        list1 = triangulation(pgon1);
        list2 = triangulation(pgon2);

        pgon = [];

        for i = 1:length(list1)
            for j = 1:length(list2)
                pgon = pgon | convHull(list1{i}, list2{j});
            end
        end

        warning(w);
    end

    function I = interval(pgon)
        % compute an interval enclosure of the polygon

        V = pgon.set.Vertices;
        I = interval(min(V, [], 1)', max(V, [], 1)');
    end

    function P = polytope(pgon)
    % convert a convex polygon to a polytope object
        if isConvex(pgon)
            P = polytope(pgon.set.Vertices');
        else
            P = polytope(convHull(pgon)');
        end

    end

    function p = randPoint(pgon, varargin)
        % draw random points within the polygon

        % parse input arguments
        [N, type] = setDefaultValues({1, 'standard'}, varargin);

        % different types of random points
        if strcmp(type, 'standard')

            p = zeros(2, N);

            list = triangulation(pgon);
            cnt = 1;
            cnt_ = 1;

            while cnt <= N

                % get random point by interpolation between the vertices
                V = list{cnt_}.set.Vertices;
                d = rand(3, 1);
                d = d ./ sum(d);
                p(:, cnt) = sum(V'*diag(d), 2);

                % update counter
                cnt = cnt + 1;
                cnt_ = cnt_ + 1;
                if cnt_ > length(list)
                    cnt_ = 1;
                end
            end

        elseif strcmp(type, 'extreme')

            % return all extreme point
            if ischar(N) && strcmp(N, 'all')

                p = pgon.set.Vertices';

            else

                V = pgon.set.Vertices;

                if N <= size(V, 1)
                    ind = randperm(size(V, 1));
                    ind = ind(1:N);
                    p = V(ind, :)';
                else
                    p = [V', zeros(2, size(V, 1)-N)];
                    for i = size(V, 1) + 1:N
                        ind = randi(size(V, 1)-1, [1, 1]);
                        d = rand(2, 1);
                        d = d / sum(d);
                        v1 = V(ind, :)';
                        v2 = V(ind+1, :)';
                        p(:, i) = d(1) * v1 + d(2) * v2;
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
            if length(convhull(pgon.set.Vertices)) - 1 == size(pgon.set.Vertices, 1)
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
        list = cell(size(T.ConnectivityList, 1), 1);

        for i = 1:size(T.ConnectivityList, 1)
            V = T.Points(T.ConnectivityList(i, :), :);
            list{i} = polygon(V(:, 1), V(:, 2));
        end
    end

    function pgon = simplify(pgon, varargin)
        % enclose the polygon by a simpler polygon with less vertices

        if nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs', 2));
        end
        % default values
        tol = setDefaultValues({0.01}, varargin);
        % check input arguments
        inputArgsCheck({{pgon, 'att', 'polygon'}, ...
            {tol, 'att', 'numeric', {'scalar', 'nonnegative'}}});

        % simplify polygon boundary using the Douglar-Peucker algorithm
        V = pgon.set.Vertices;
        V_ = douglasPeucker(V', tol);
        pgon = polygon(V_(1, :), V_(2, :));

        % enlarge the simplified polygon so that it is guaranteed that it
        % encloses the original polygon
        temp = polygon([-tol, -tol, tol, tol], [-tol, tol, tol, -tol]);
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
                for j = i + 1:length(list)

                    tmp = list{i} | list{j};

                    if isConvex(tmp)
                        list{i} = tmp;
                        list{j} = [];
                        stop = true;
                        finished = false;
                        break;
                    end
                end

                if stop
                    break;
                end
            end

            % remove empty list entries
            list = list(~cellfun('isempty', list));

            if finished
                break;
            end
        end
    end

    function res = isequal(pgon1, pgon2, varargin)
        % compares whether two polygons represent the same shape
        % notes: - (almost) collinear points are removed by the constructor
        %        - vertices have to be in order, but start of list can vary

        % set tolerance
        tol = setDefaultValues({1e-8}, varargin);

        % the constructor already removes collinear points, but within its
        % own tolerance -> look for collinear points up to given tolerance
        pgon1 = removeCollinearPoints(pgon1, tol);
        pgon2 = removeCollinearPoints(pgon2, tol);

        % read out vertices
        V1 = pgon1.set.Vertices;
        V2 = pgon2.set.Vertices;

        % check empty case
        if isempty(V1) && isempty(V2)
            res = true; return
        end

        % since the constructor removes collinear points, the number of
        % vertices has to be equal
        if size(V1, 1) ~= size(V2, 1)
            res = false; return
        end

        % to deal with the different start vertex, first find one vertex
        % that is part of both lists
        % -> check all occurrences of the first vertex of P1 in P2
        idx = find(all(withinTol(V1(1, :), V2, tol), 2));

        % any matching vertex found?
        if isempty(idx)
            res = false; return
        end

        % loop over all potential start vertices
        for i = 1:length(idx)
            % re-order vertices of P2
            V2_ = [V2(idx:end, :); V2(1:idx-1, :)];
            % compare in order
            if compareMatrices(V1, V2_, tol, 'equal', true)
                res = true; return
            end
        end

        % no re-ordering successful
        res = false;

    end

    function res = eq(pgon1, pgon2, varargin)
        % overloads '==' operator
        res = isequal(pgon1, pgon2, varargin{:});
    end

    function res = ne(pgon1, pgon2, varargin)
        % overloads '~=' operator
        res = ~isequal(pgon1, pgon2, varargin{:});
    end

    function pgon = removeCollinearPoints(pgon, varargin)
        % removes collinear points up to a user-defined tolerance by
        % checking the rank of two vectors computed from two subsequent
        % pairs of vertices in the list

        % set default tolerance
        tol = setDefaultValues({0}, varargin);

        % read out ordered vertices
        V = pgon.set.Vertices;

        % remove collinear vertices
        pgon = polygon(removeCollinearVertices2D(V',tol));

    end

    function res = isemptyobject(pgon)
        res = isempty(pgon.set.Vertices);
    end

    function [res,S] = representsa_(pgon,type,tol,varargin)

        switch type
            case 'emptySet'
                res = isempty(pgon.set.Vertices);
                S = [];
            otherwise
                throw(CORAerror('CORA:notSupported',...
                    'only type = ''emptySet'' supported'));

        end

    end
end

methods (Static = true)

    function pgon = generateRandom()
        % generate random polygon

        points = 10 * rand(1) * (-1 + 2 * rand(2, 100)) + 10 * (-1 + 2 * rand(2, 1));
        pgon = polygon.enclosePoints(points);
    end

    function pgon = enclosePoints(points)
        % enclose point cloud with polygon

        ind = boundary(points', 0.5);
        pgon = polygon(points(1, ind), points(2, ind));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
