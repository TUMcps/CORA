classdef conZonotope < contSet
% conZonotope - object constructor for constrained zonotopes [1]
%
% Description:
%    This class represents constrained zonotope objects defined as
%    {c + G * beta | ||beta||_Inf <= 1, A * beta <= b}.
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
%    obj - generated conZonotope object
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
% Last update:   19-March-2021 (MW, errConstructor)
% Last revision: 02-May-2020 (MW, methods list, rewrite methods(hidden),
%                                 add property validation)

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    % center and generators  x = Z(:,1) + Z(:,2:end)*ksi; |ksi| <= 1
    % format:       matrix
    Z (:,:) {mustBeNumeric,mustBeFinite} = [];
    
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
        
        Z = []; A = []; b = [];
        
        if nargin == 1
            % copy constructor
            if isa(varargin{1},'conZonotope')
                obj = varargin{1};
            else
                Z = varargin{1};
            end
        
        elseif nargin == 2
            if ~isvector(varargin{1}) || length(varargin{1}) ~= size(varargin{2},1)
                [id,msg] = errConstructor(); error(id,msg);
            end
            Z = [varargin{1},varargin{2}];
            
        elseif nargin == 3
            Z = varargin{1};
            A = varargin{2};
            b = varargin{3};
            
        elseif nargin == 4
            if ~isvector(varargin{1}) || length(varargin{1}) ~= size(varargin{2},1)
                [id,msg] = errConstructor(); error(id,msg);
            end
            Z = [varargin{1},varargin{2}];
            A = varargin{3};
            b = varargin{4};
            
        elseif nargin > 4
            % too many input arguments
            [id,msg] = errConstructor('Too many input arguments.'); error(id,msg);
        end
        
        if ~isempty(A) && ~isempty(b)
            % check correctness of A and b and w.r.t Z
            if ~isvector(b) ... % b is a vector ?
                    || size(A,2) ~= (size(Z,2)-1) ... % A fits Z (gens)?
                    || size(A,1) ~= length(b) % A fits b ?
                [id,msg] = errConstructor(); error(id,msg);
            end
        end
        
        % store object properties
        obj.Z = Z; obj.A = A; obj.b = b;
        
        % set parent object properties
        obj.dimension = size(obj.Z,1);
    end
    
    
    % methods in seperate files
    res = and(obj,S);
    cZ = cartProd(cZ1,cZ2);
    res = center(obj);
    cZ = convHull(cZ1,varargin);
    cPZ = conPolyZono(cZ);
    d = dim(obj);
    res = deleteZeros(obj);
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
    res = reduceConstraints(obj,varargin);
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

end

%------------- END OF CODE -------