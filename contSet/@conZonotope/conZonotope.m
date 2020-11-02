classdef conZonotope < zonotope
% conZonotope - object constructor for constrained zonotopes [1]
%
% Syntax:
%       obj = conZonotope(c,G)
%       obj = conZonotope(c,G,A,b)
%       obj = conZonotope(Z)
%       obj = conZonotope(Z,A,b)
%
% Inputs:
%    c - center vector of the zonotope
%    G - generator matrix of the zonotope
%    Z - matrix containing zonotope center and generators Z = [c,G]
%    A - constraint matrix A*ksi = b
%    b - constraint vector A*ksi = b
%
% Outputs:
%    obj - Generated Object
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono = conZonotope(Z,A,b);
%    plotZono(cZono);
%
% References:
%    [1] Scott, Joseph K., et al. "Constrained zonotopes:
%           A new tool for set-based estimation and fault detection."
%           Automatica 69 (2016): 126-136.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval
%
% Author:        Dmitry Grebenyuk, Mark Wetzlinger
% Written:       03-September-2017
% Last update:   ---
% Last revision: 02-May-2020 (MW, methods list, rewrite methods(hidden),
%                                 add property validation)

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    % constraint A*ksi = b; |ksi| <= 1
    % format:       matrix
    A (:,:) {mustBeNumeric,mustBeFinite} = [];
    
    % format:       column vector
    b (:,1) {mustBeNumeric,mustBeFinite} = [];
    
    % the value of ksi at vertexes
    % format:       column vector
    ksi (:,:) {mustBeNumeric,mustBeFinite} = [];
    
    % R = [rho_l, rho_h] (A.3)
    % format:       column vector
    R = [];
    
end
    
methods
    
    % class constructor
    function obj = conZonotope(varargin)
        
        A = [];
        b = [];
        
        if nargin == 0
            Z = [];
        elseif nargin == 1
            Z = varargin{1};
        elseif nargin == 2
            Z = [varargin{1},varargin{2}];
        elseif nargin == 3
            Z = varargin{1};
            A = varargin{2};
            b = varargin{3};
        elseif nargin == 4
            Z = [varargin{1},varargin{2}];
            A = varargin{3};
            b = varargin{4};
        else
            error('This class takes at max 4 inputs.')
        end
        
        % create a zonotope
        obj@zonotope(Z);
        
        % check if A and b fit the zonotope
        % note: this check is necessary, otherwise run-time errors will occur
        if ~isempty(A) && ~isempty(b)
            if size(A,2) ~= (size(obj.Z,2)-1)
                error("A has to be of proper dimension.");
            end
            if length(b) ~= size(A,1)
                error("b has to be of proper dimension.");
            end
        end
        
        % store the matrices for the constraints
        obj.A = A;
        obj.b = b;
    end
    
    
    % methods in seperate files
    res = and(obj,S);
    cZ = cartProd(cZ1,cZ2);
    res = center(obj);
    cZ = convHull(cZ1,varargin);
    d = dim(obj);
    display(obj);
    cZ = enclose(varargin);
    res = in(obj1,obj2,varargin);
    res = interval(obj);
    res = intervalMultiplication(obj,I);
    res = isempty(obj);
    res = isFullDim(obj);
    res = isIntersecting(obj1,obj2,varargin);
    res = mptPolytope(obj);
    cZ = mtimes(factor1,factor2);
    cZ = or(cZ1, varargin);
    handle = plot(obj,varargin);
    handle = plotZono(obj,varargin);
    cZ = plus(summand1,summand2);
    pZ = polyZonotope(obj);
    obj = project(obj,dims);
    cZquad = quadMap(varargin);
    res = reduce(obj,method,orderG,varargin);
    res = rescale(obj,varargin);
    cZsplit = split(obj,varargin);
    [val,x,ksi] = supportFunc(obj,dir,varargin);
    V = vertices(obj);
    res = zonoBundle(obj);
    res = zonotope(obj,varargin);    
             
end

methods (Static = true)
    cZ = generateRandom(varargin) % generate random constrained zonotope
end

% prevent unintentional usage of superclass methods
methods(Hidden)
    
    function abs(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function box(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function capsule(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function constrSat(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function containsPoint(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function deleteAligned(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function deleteZeros(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function dominantDirections(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function ellipsoid(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function encloseMany(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function enclosingPolytope(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function enlarge(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function eventFcn(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function exactPlus(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function filterOut(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function generatorLength(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function generators(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function get(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function halfspace(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function intersectZonoStrip(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function isequal(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function isInterval(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function minnorm(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function minus(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function norm(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function normbound(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function orthVectors(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function polygon(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function polytope(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function quadMap_parallel(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function radius(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function rank(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function reduceUnderApprox(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function rotate(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function sampleBox(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function set(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function splitFirstGen(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function taylm(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function tensorMultiplication(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function tensorMultiplication_zono(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function underapproximate(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function volume(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
    function volumeRatio(varargin)
        error('This operation is not implemented for class "conZonotope"!'); 
    end
    
end


end

%------------- END OF CODE -------