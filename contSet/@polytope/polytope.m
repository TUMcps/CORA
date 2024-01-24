classdef polytope < contSet
% polytope - object constructor for polytope objects
%
% Description:
%    This class represents polytope objects defined as
%      { x | A*x <= b }. 
%    For convenience, equality constraints
%      { x | A*x <= b, Ae*x == be}
%    can be added, too.
%    Note: A polytope without any constraints represents R^n.
%    Note: A polytope instantiated without input arguments is the empty set.
%
% Syntax:
%    P = polytope(V)
%    P = polytope(A,b)
%    P = polytope(A,b,Ae,be)
%
% Inputs:
%    V - (n x p) array of vertices (gets converted to halfspace representation)
%    A - (n x m) matrix for the inequality representation
%    b - (n x 1) vector for the inequality representation
%    Ae - (k x l) matrix for the equality representation
%    be - (k x 1) vector for the equality representation
%
% Outputs:
%    obj - generated polytope object
%
% Example: 
%    A = [1 0 -1 0 1; 0 1 0 -1 1]';
%    b = [3; 2; 3; 2; 1];
%    poly = polytope(A,b);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Authors:       Viktor Kotsev, Mark Wetzlinger, Tobias Ladner
% Written:       25-April-2022
% Last update:   01-December-2022 (MW, add CORAerrors, checks)
%                12-June-2023 (MW, add hidden properties)
%                08-December-2023 (MW, handle -Inf/Inf offsets)
%                01-January-2024 (MW, different meaning of fully empty obj)
% Last revision: 25-July-2023 (MW, restructure constructor)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)

    % Inequality description { x | A*x <= b }
    A = [];
    b = [];

    % Affine set description { x | Aeq*x == be }
    Ae = [];
    be = [];
    
end

properties (SetAccess = protected, GetAccess = public)
    % Vertex description
    V;

    % emptiness
    emptySet;

    % full-dimensionality
    fullDim;

    % boundedness
    bounded;

    % minimal halfspace representation
    minHRep;

    % minimal vertex representation
    minVRep;

end

methods
    function obj = polytope(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 0. init hidden properties
        obj.emptySet = setproperty();
        obj.fullDim = setproperty();
        obj.bounded = setproperty();
        obj.minHRep = setproperty();
        obj.minVRep = setproperty();
        obj.V = setproperty();
    
        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'polytope')
            % read out polytope
            P = varargin{1};
            % copy properties
            obj.A = P.A; obj.b = P.b;
            obj.Ae = P.Ae; obj.be = P.be;
            % copy hidden properties
            obj.emptySet.val = P.emptySet.val;
            obj.fullDim.val = P.fullDim.val;
            obj.bounded.val = P.bounded.val;
            obj.minHRep.val = P.minHRep.val;
            obj.minVRep.val = P.minVRep.val;
            obj.V.val = P.V.val;
            return
        end

        % 2. parse input arguments: varargin -> vars
        [A,b,Ae,be,V] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(A,b,Ae,be,V,nargin);

        % 4. compute properties and hidden properties
        [A,b,Ae,be,V] = aux_computeProperties(A,b,Ae,be,V);
        [empty,bounded,fullDim,minHRep,minVRep,V] = ...
            aux_computeHiddenProperties(A,b,Ae,be,V);

        % 4a. assign properties
        obj.A = A;
        obj.b = b;
        obj.Ae = Ae;
        obj.be = be;
        obj.V.val = V;
        obj.emptySet.val = empty;
        obj.bounded.val = bounded;
        obj.fullDim.val = fullDim;
        obj.minHRep.val = minHRep;
        obj.minVRep.val = minVRep;
        
    end

end

methods (Static = true)
    P = generateRandom(varargin) % generate random polytope
    P = enclosePoints(points,varargin) % enclose point cloud with polytope
    P = empty(n) % instantiate empty polytope
    P = Inf(n) % instantiate polytope representing R^n
end

end


% Auxiliary functions -----------------------------------------------------

function [A,b,Ae,be,V] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4));
    elseif nargin == 3
        % no syntax for three input arguments
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Constructor for class polytope requires 0, 1, 2, or 4 input arguments.'));
    end

    % no input arguments
    if nargin == 0
        A = []; b = []; Ae = []; be = []; V = [];
        return
    end

    % read out arguments
    if nargin == 1
        % vertices as input argument
        V = varargin{1};
        A = []; b = []; Ae = []; be = [];
    else
        % halfspaces as input arguments
        [A,b,Ae,be] = setDefaultValues({[],[],[],[]},varargin);
        V = [];
    end

end

function aux_checkInputArgs(A,b,Ae,be,V,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % check numeric type of V
        if ~isnumeric(V)
            if any(any(isnan(V)))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Vertices have to be non-nan.'));
            elseif size(V,1) > 1 && any(any(isinf(V)))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'nD vertices for n > 1 have to be finite.'));
            end
        end

        % check b, be (get later reshaped to column vector)
        if ~isempty(b) && ~isvector(b)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Arguments "b" has to be column vector.'));
        end
        if ~isempty(be) && ~isvector(be)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Arguments "be" has to be column vector.'));
        end

        % check number of inequality constraints
        if ~isempty(A) && size(A,1) ~= length(b)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Number of rows does not hold between arguments "A", "b".'));
        end
        
        % check for empty argument
        if isempty(A) && ~isempty(b) && size(A,1) ~= length(b)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Number of rows does not hold between arguments "A", "b".'));
        end

        % check number of equality constraints
        if ~isempty(Ae) && size(Ae,1) ~= length(be)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Number of rows does not hold between arguments "Ae", "be".'));
        end

        % same dimension if both equality and inequality constraints given
        if ~isempty(A) && ~isempty(Ae) && size(A,2) ~= size(Ae,2)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Number of columns does not hold between arguments "A", "Ae".'));
        end
    end

end

function [A,b,Ae,be,V] = aux_computeProperties(A,b,Ae,be,V)

    % dimension
    n = max([size(A,2),size(Ae,2),size(V,1)]);

    % offsets must be column vectors
    b = reshape(b,[],1);
    be = reshape(be,[],1);

    % compute halfspace representation
    if ~isempty(V)
        % for n=1, remove redundancies
        if n == 1
            V = [min(V), max(V)];
            if withinTol(V(1),V(2),eps)
                V = V(1);
            end
        end
        [A,b,Ae,be] = computeHRep(V);
    else
        V = [];
    end

    % empty constraint matrices should have correct dimension (necessary
    % for matrix concatenation)
    if isempty(A)
        A = zeros(0,n);
    end
    if isempty(Ae)
        Ae = zeros(0,n);
    end

    % remove inequality constraints with Inf in offset (trivially fulfilled)
    idxRemove = isinf(b) & sign(b) == 1;
    A(idxRemove,:) = []; b(idxRemove) = [];

end

function [empty,bounded,fullDim,minHRep,minVRep,V] = ...
            aux_computeHiddenProperties(A,b,Ae,be,V)

    % init hidden properties as unknown
    empty = []; bounded = []; fullDim = []; minHRep = []; minVRep = [];    

    % dimension
    n = max([size(A,2),size(Ae,2),size(V,1)]);

    % representations for 1D polytopes always minimal
    if n == 1 && ~isempty(V)
        minHRep = true;
    end

    % check if instantiated via vertices
    if ~isempty(V)
        % cannot be empty
        empty = false;

        % only 1 vertex -> minimal V-representation
        minVRep = size(V,2) == 1;

        % check if 1D
        if n == 1
            % Inf values for vertices only supported for 1D
            if any(any(isinf(V)))
                bounded = false;
                fullDim = true;
            else
                bounded = true;
                fullDim = size(V,2) > 1;
            end
            
        else
            % nD -> has to bounded
            bounded = true;
            % easy checks for degeneracy
            if size(V,2) <= n
                % full-dimensionality requires at least n+1 vertices
                fullDim = false;
            else
                [~,S,~] = svd(V);
                fullDim = size(V,1) == nnz(~withinTol(S,0,1e-12));
            end
        end

    elseif isempty(A) && isempty(Ae)
        % no constraints
        empty = false;
        bounded = false;
        fullDim = true;
        minHRep = true;
        % do not compute -Inf/Inf vertices here...
        V = [];
        minVRep = [];

    else
        % equality constraint with -Inf or Inf in offset OR inequality
        % constraint with -Inf in offset -> empty polytope
        if any(isinf(be)) || ( any(isinf(b) & sign(b) == -1) )
            empty = true;
            bounded = true;
            fullDim = false;
            % only a single infeasible constraint required to represent an
            % empty set
            minHRep = length(be) + length(b) == 1;
            % init no vertices (which is the minimal representation)
            V = zeros(n,1);
            minVRep = true;
        end

    end

end

% ------------------------------ END OF CODE ------------------------------
