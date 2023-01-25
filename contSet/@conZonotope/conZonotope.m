classdef conZonotope < contSet
% conZonotope - object constructor for constrained zonotopes [1]
%
% Description:
%    This class represents constrained zonotope objects defined as
%    {c + G * beta | ||beta||_Inf <= 1, A * beta = b}.
%
% Syntax:
%    obj = conZonotope(c,G)
%    obj = conZonotope(c,G,A,b)
%    obj = conZonotope(Z)
%    obj = conZonotope(Z,A,b)
%
% Inputs:
%    c - center vector of the zonotope
%    G - generator matrix of the zonotope
%    Z - matrix containing zonotope center and generators Z = [c,G]
%    A - constraint matrix A*beta = b
%    b - constraint vector A*beta = b
%
% Outputs:
%    obj - generated conZonotope object
%
% Example: 
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%    plot(cZ);
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

% Author:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:      03-September-2017
% Last update:  19-March-2021 (MW, error messages)
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:02-May-2020 (MW, methods list, rewrite methods(hidden),
%                                 add property validation)

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    % center and generators  x = Z(:,1) + Z(:,2:end)*beta; |beta| <= 1
    % format:       matrix
    Z = [];
    
    % constraint A*beta = b; |beta| <= 1
    % format:       matrix
    A = [];
    
    % format:       column vector
    b = [];
    
    % the value of beta at vertexes
    % format:       column vector
    ksi (:,:) {mustBeNumeric,mustBeFinite} = [];
    
    % R = [rho_l, rho_h] (A.3)
    % format:       column vector
    R = [];
    
end
    
methods
    
    % class constructor
    function obj = conZonotope(varargin)

        % parse input
        if nargin > 4
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end
        
        Z = []; A = []; b = [];
        inputChecks = {};
        
        if nargin == 1 && isa(varargin{1},'conZonotope')
                % copy constructor
                obj = varargin{1};
                return;
            
        elseif nargin == 1 || nargin == 3
            Z = varargin{1};
            varargin = varargin(2:end);

            inputChecks = {{Z, 'att', 'numeric', 'matrix'}};
        
        elseif nargin == 2 || nargin == 4
            c = varargin{1};
            G = varargin{2};
            varargin = varargin(3:end);

            inputChecks = { ...
                {c, 'att', 'numeric', {'finite', 'column'}};
                {G, 'att', 'numeric', {'finite', 'matrix'}};
            };
            inputArgsCheck(inputChecks);

            % check dimensions
            if length(c) ~= size(G,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimensions of the center and the generator matrix do not match.')); 
            end
            Z = [c,G];
        end

        [A, b] = setDefaultValues({[], []}, varargin);

        % check dimensions
        inputArgsCheck([ ...
            inputChecks; ...
            {{A, 'att', 'numeric', {'finite', 'matrix'}}};
            {{b, 'att', 'numeric', {'finite'}}};
        ])
        
        if ~isempty(A) && ~isempty(b)
            % check correctness of A and b and w.r.t Z
            if ~isvector(b) % b is a vector ?
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The constraint offset has to be a vector.'));
            elseif size(A,2) ~= (size(Z,2)-1) % A fits Z (gens)?
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimensions of the generator matrix and the constraint matrix do not match.'));
            elseif size(A,1) ~= length(b) % A fits b ?
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimensions of the constraint matrix and the constraint offset do not match.'));
            end
        end
        
        % store object properties
        obj.Z = Z; obj.A = A; obj.b = b;
        
        % set parent object properties
        obj.dimension = size(obj.Z,1);
    end
    
    
    % methods in seperate files
    res = and(cZ,S)
    cZ = cartProd(cZ,S)
    c = center(cZ)
    res = conIntersect(cZ1,cZ2,M)
    cPZ = conPolyZono(cZ)
    res = contains(cZ,S,varargin)
    cZ = convHull(cZ,varargin)
    res = cubMap(cZ,varargin)
    cZ = deleteZeros(cZ)
    n = dim(cZ)
    display(cZ);
    cZ = enclose(cZ,varargin)
    G = generators(cZ)
    cZ = intersectStrip(cZ,C,phi,y,varargin)
    I = interval(cZ)
    cZ = intervalMultiplication(cZ,I)
    res = isempty(cZ)
    res = isFullDim(cZ)
    res = isIntersecting(cZ,S,varargin)
    cZ = minkDiff(cZ1,S,varargin)
    P = mptPolytope(cZ)
    cZ = mtimes(factor1,factor2)
    cZ = or(cZ1,varargin)
    han = plot(cZ,varargin)
    han = plotZono(cZ,varargin)
    cZ = plus(summand1,summand2)
    pZ = polyZonotope(cZ)
    cZ = project(cZ,projDim)
    cZ = quadMap(cZ,varargin)
    p = randPoint(cZ,varargin)
    cZ = reduce(cZ,method,order,varargin)
    cZ = reduceConstraints(cZ,varargin)
    cZ = rescale(cZ,varargin)
    cZsplit = split(cZ,varargin)
    [val,x,ksi] = supportFunc(cZ,dir,varargin)
    V = vertices(cZ)
    zB = zonoBundle(cZ)
    Z = zonotope(cZ,varargin)
             
end

methods (Static = true)
    cZ = generateRandom(varargin) % generate random constrained zonotope
end

end

%------------- END OF CODE --------------