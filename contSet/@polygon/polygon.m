classdef polygon < contSet
% polygon class - vertices of a nonconvex set in 2D
%    based on Matlab's polyshape class
%
% Syntax:
%    obj = polygon(x,y,varargin)
%    obj = polygon(V,varargin)
%    obj = polygon(set,varargin)
%
% Inputs:
%    x - vector with x coordinates of the polygon vertices
%    y - vector with y coordinates of the polygon vertices
%    V - vertices (2-dimensional)
%    set - polyshape object
%    varargin - name-value pairs for polyshape constructor
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
% See also: polyshape

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       13-March-2020
% Last update:   09-May-2023 (TL, constructor clean up)
%                13-March-2024 (TL, enable polyshape NVpairs in constructor)
%                08-October-2024 (TL, made proper class within contSet)
%                23-October-2024 (TL, added V property)
%                13-February-2024 (TL, added nrOfRegions/nrOfHoles properties)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    set = []
    % dummy properties
    x = []; y = [];
end

properties (Access=private)
    V;          % vertices (if present)
    TOL = 1e-8; % for buffered degenerate polygons
end

properties (Dependent)
    nrOfRegions; nrOfHoles;
end

methods

    % class constructor
    function obj = polygon(varargin)

        % parse input
        [set, x, y, NVpairs] = aux_parseInputArgs(obj, varargin{:});

        % check input
        aux_checkInputArgs(set,x,y,NVpairs)

        % postprocessing
        obj = aux_postprocessing(obj,set,x,y,NVpairs);

        % 6. set precedence (fixed)
        obj.precedence = 15;
    end

end

methods (Static = true)
    % static functions
    pgon = generateRandom(varargin);
    pgon = enclosePoints(points);
    pgon = empty(n)
end

methods (Access = protected)
    % protected functions
    [abbrev,printOrder] = getPrintSetInfo(S)
end


% getter & setter ---------------------------------------------------------

methods
    function x = get.x(pgon)
        % get x
        V = vertices_(pgon);
        x = V(1,:);
    end

    function y = get.y(pgon)
        % get y
        V = vertices_(pgon);
        y = V(2,:);
    end

    function nrOfRegions = get.nrOfRegions(pgon)
        % get number of regions
        nrOfRegions = pgon.set.NumRegions;
    end

    function nrOfHoles = get.nrOfHoles(pgon)
        % get number of holes
        nrOfHoles = pgon.set.NumHoles;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [set, x, y, NVpairs] = aux_parseInputArgs(obj, varargin)
    % default values
    set = []; x = []; y = []; NVpairs = {};

    % empty constructor
    if nargin == 1
        set = polyshape();
        return
    end
    
    % polygon(set)
    if nargin == 2
        set = varargin{1};
        return
    end

    % nargin > 2 ---

    % check if set or (x,y) pair is given
    x = varargin{1};
    y = varargin{2};

    % check if second parameter is numeric or first name-value pair
    if isnumeric(y)
        % (x,y) given
        NVpairs = varargin(3:end);
    else
        % 'y' is first name of first name-value pair
        set = x;
        x = []; y = [];
        NVpairs = varargin(2:end);
    end
end

function aux_checkInputArgs(set,x,y,NVpairs)
    % check set/V
    inputArgsCheck({{set, 'att', {'polygon', 'polyshape', 'numeric'}}})
    if isnumeric(set) 
        inputArgsCheck({{set, 'att', 'numeric'}})
        if size(set,1) > 0 && size(set,1) ~= 2
            % 'set' is a matrix containing vertices
            throw(CORAerror("CORA:wrongValue",'first','Given vertices should be two dimensional.'))
        end
    end

    % check (x,y)
    inputArgsCheck({{x, 'att', 'numeric'}, {y, 'att', 'numeric'}})
    if ~isempty(x) && (~all(size(x) == size(y)) || ~isvector(x))
        throw(CORAerror("CORA:wrongInputInConstructor", ...
            'Given vertices x,y need to be vectors of the same length.'))
    end

    % check if finite (but nan is ok)
    if ~all(~isinf(x)) || ~all(~isinf(y)) || (isnumeric(set) && ~all(~isinf(set),"all"))
        throw(CORAerror("CORA:wrongInputInConstructor", ...
            'Given vertices cannot be infinite.'))
    end

end

function obj = aux_postprocessing(obj,set,x,y,NVpairs)

    % polyshape given?
    if isa(set,'polyshape')
        obj.set = set;
        return
    end

    % polygon given?
    if isa(set,'polygon')
        obj.set = set.set;
        return
    end

    % 'set' is a matrix
    if isnumeric(set) && ~isempty(set)
        % rewrite to x and y
        x = set(1,:); y = set(2,:);
    end

    % empty x,y?
    if isempty(x)
        obj.set = polyshape();
        return
    end

    % x,y given
    x = reshape(x,1,[]);
    y = reshape(y,1,[]);

    % init polyshape
    w = warning; warning off;
    set = polyshape(x,y,NVpairs{:});
    warning(w);

    % read vertices
    V = set.Vertices';

    % check if vertices are given
    tol = obj.TOL;
    while isempty(V)
        % make proper set by adding some buffer
        tol = tol * 10;
        tol_x = [tol,tol,-tol,-tol];
        tol_y = [tol,-tol,-tol,tol];

        set = polyshape();
        for i=1:numel(x)
            % buffered vertex
            v = polyshape(x(i)+tol_x,y(i)+tol_y);

            % compute convex hull of set and buffered vertex
            set = convhull(union(set, v));

            % should contain each subset
            if ~overlaps(set,v)
                % restart with higher tol
                set = polyshape();
                break
            end
        end
        
        % remove holes
        set = rmholes(set);
        V = [x;y];
    end
    
    % set properties
    obj.set = set;
    obj.TOL = tol;
    obj.V = V;
end

% ------------------------------ END OF CODE ------------------------------
