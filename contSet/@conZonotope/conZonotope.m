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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger, Tobias Ladner
% Written:       03-September-2017
% Last update:   19-March-2021 (MW, error messages)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                29-March-2023 (TL, optimized constructor)
%                13-September-2023 (TL, replaced Z property with c and G)
% Last revision: 02-May-2020 (MW, methods list, rewrite methods(hidden), add property validation)
%                16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    % center and generators  x = c + G*beta; |beta| <= 1
    c, G;
    Z = []; % legacy Z = [c,g_1,...,g_p]
    
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

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'conZonotope')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [c,G,A,b] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(c,G,A,b,nargin);

        % 4. compute properties
        [c,G,A,b] = aux_computeProperties(c,G,A,b);

        % 5. assign properties
        obj.c = c;
        obj.G = G;
        obj.A = A;
        obj.b = b;

    end
    
    
    % methods in seperate files
    c = center(cZ)
    res = conIntersect(cZ1,cZ2,M)
    cPZ = conPolyZono(cZ)
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
    res = isFullDim(cZ)
    cZ = minkDiff(cZ1,S,varargin)
    P = polytope(cZ)
    cZ = mtimes(factor1,factor2)
    cZ = or(cZ1,varargin)
    han = plot(cZ,varargin)
    han = plotZono(cZ,varargin)
    cZ = plus(summand1,summand2)
    pZ = polyZonotope(cZ)
    cZ = project(cZ,projDim)
    cZ = quadMap(cZ,varargin)
    cZ = reduce(cZ,method,order,varargin)
    cZ = reduceConstraints(cZ,varargin)
    cZ = rescale(cZ,varargin)
    [res,S] = representsa_(cZ,type,tol,varargin)
    cZsplit = split(cZ,varargin)
    zB = zonoBundle(cZ)
    Z = zonotope(cZ,varargin)             
end

methods (Static = true)
    cZ = generateRandom(varargin) % generate random constrained zonotope
    cz = empty(n) % instantiates an empty constrained zonotope
end


% getter & setter ---------------------------------------------------------

methods
    function obj = set.G(obj,G)
        % fix dimension if empty
        if isempty(G)
            G = zeros(dim(obj),0);
        end
        obj.G = G;
    end

    % getter & setter for legacy Z property
    function Z = get.Z(obj)
        warning(['CORA: The property conZonotope.Z is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conZonotope.c and conZonotope.G instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        Z = [obj.c, obj.G];
    end

    function obj = set.Z(obj, Z)
        warning(['CORA: The property conZonotope.Z is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conZonotope.c and conZonotope.G instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        if ~isempty(Z)
            obj.c = Z(:,1);
            obj.G = Z(:,2:end);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G,A,b] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4));
    end

    % set default values depending on nargin
    if nargin == 1 || nargin == 3
        % only center given, or [c,G] with A and b
        [c,A,b] = setDefaultValues({[],[],[]},varargin);
        if size(varargin{1},2) > 0
            G = c(:,2:end);
            c = c(:,1);
        else
            G = [];
        end
    elseif nargin == 2 || nargin == 4
        % c,G or c,G,A,b given
        [c,G,A,b] = setDefaultValues({[],[],[],[]},varargin);
    end

end

function aux_checkInputArgs(c,G,A,b,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        if n_in == 1 || n_in == 3
            inputChecks = { ...
                {c, 'att', 'numeric', {'finite'}};
                {G, 'att', 'numeric', {'finite', 'matrix'}};
            };

        elseif n_in == 2 || n_in == 4
            % check whether c and G fit together to avoid bad message
            if ~isempty(c) && ~isvector(c)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The center has to be a column vector.')); 
            elseif ~isempty(G) && length(c) ~= size(G,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimensions of the center and the generator matrix do not match.')); 
            end
            inputChecks = {{[c,G], 'att', 'numeric', 'finite'}};
        end

        inputArgsCheck([inputChecks; ...
            {{A, 'att', 'numeric', {'finite', 'matrix'}};
            {b, 'att', 'numeric', 'finite'}}])

        % check correctness of A and b, also w.r.t G
        if ~isempty(A) && ~isempty(b)
            if ~isvector(b) % b is a vector?
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The constraint offset has to be a vector.'));
            elseif size(A,2) ~= (size(G,2)) % A fits G?
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimensions of the generator matrix and the constraint matrix do not match.'));
            elseif size(A,1) ~= length(b) % A fits b?
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The dimensions of the constraint matrix and the constraint offset do not match.'));
            end
        end
        
    end

end

function [c,G,A,b] = aux_computeProperties(c,G,A,b)

    % if G is empty, set correct dimension
    if isempty(G)
        G = zeros(size(c,1),0);
    end

    % if no constraints, set correct dimension
    if isempty(A)
        A = zeros(0,size(G,2));
        b = zeros(0,1);
    end

end

% ------------------------------ END OF CODE ------------------------------
