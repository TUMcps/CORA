classdef polygon < contSet
% polygon class - vertices of a nonconvex set in 2D
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    set = []
    % dummy properties
    x = []; y = [];
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
    pgon = generateRandom(varargin);
    pgon = enclosePoints(points);
    pgon = empty(n)
end

methods (Access = protected)
    [abbrev,printOrder] = getPrintSetInfo(S)
end


% getter & setter ---------------------------------------------------------

methods
    function x = get.x(pgon)
        V = vertices(pgon);
        x = V(1,:);
    end

    function y = get.y(pgon)
        V = vertices(pgon);
        y = V(2,:);
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
    w = warning; warning off;
    set = polyshape(x,y,NVpairs{:});
    warning(w);

    % check if vertices are given
    V = set.Vertices;
    if isempty(V)
        % make proper set by adding some buffer
        tol = 1e-7;
        tol_x = [tol,tol,-tol,-tol];
        tol_y = [tol,-tol,-tol,tol];
        
        % compute boundary
        x = [reshape(x,1,[]) x(1)];
        y = [reshape(y,1,[]) y(1)];

        set = polyshape();
        for i=1:(numel(x)-1)
            p1 = polyshape(x(i)+tol_x,y(i)+tol_y);
            p2 = polyshape(x(i+1)+tol_x,y(i+1)+tol_y);

            set = union(set, convhull(union(p1,p2)));
        end
        
        % remove holes
        set = rmholes(set);
    end
    
    obj.set = set;
    return
end

% ------------------------------ END OF CODE ------------------------------
