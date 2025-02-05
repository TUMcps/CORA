classdef polytope < contSet
% polytope - object constructor for polytope objects
%
% Description:
%    This class represents polytope objects defined as (halfspace
%    representation)
%      { x | A*x <= b }. 
%    For convenience, equality constraints
%      { x | A*x <= b, Ae*x == be }
%    can be added, too.
%    Alternatively, polytopes can be defined as (vertex representation)
%      { sum_i a_i v_i | sum_i a_i = 1, a_i >= 0 }
%    Note: A polytope without any constraints represents R^n.
%    Note: A polytope instantiated without input arguments is the empty set.
%
% Syntax:
%    P = polytope(V)
%    P = polytope(A,b)
%    P = polytope(A,b,Ae,be)
%
% Inputs:
%    V - (n x p) array of vertices
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
%    P = polytope(A,b);
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
%                13-March-2024 (TL, check if input is numeric)
%                16-July-2024 (MW, allow separate usage of VRep/HRep)
% Last revision: 25-July-2023 (MW, restructure constructor)

% ------------------------------ BEGIN CODE -------------------------------

% halfspace representation and vertex representation
properties (SetAccess = private, GetAccess = private)

    % Inequality description { x | A*x <= b }
    A_ = [];
    b_ = [];

    % Affine set description { x | Ae*x == be }
    Ae_ = [];
    be_ = [];

    % Vertex description { sum_i a_i v_i | sum_i a_i = 1, a_i >= 0 }
    V_;
    
end

properties (Dependent)
    
    % public access for constraints and vertices, see above
    A;
    b;
    Ae;
    be;
    V;

end

% set properties
properties (SetAccess = protected, GetAccess = public)

    % halfspace representation given?
    isHRep;

    % vertex representation given?
    isVRep;

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
        assertNarginConstructor([1,2,4],nargin);

        % 0. init setproperty properties
        obj.A_ = setproperty();
        obj.b_ = setproperty();
        obj.Ae_ = setproperty();
        obj.be_ = setproperty();
        obj.V_ = setproperty();

        obj.isHRep = setproperty();
        obj.isVRep = setproperty();
        obj.emptySet = setproperty();
        obj.fullDim = setproperty();
        obj.bounded = setproperty();
        obj.minHRep = setproperty();
        obj.minVRep = setproperty();
    
        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'polytope')
            % read out polytope
            P = varargin{1};

            % copy properties
            obj.A_.val = P.A_.val; obj.b_.val = P.b_.val;
            obj.Ae_.val = P.Ae_.val; obj.be_.val = P.be_.val;
            obj.V_.val = P.V_.val;
            obj.precedence = P.precedence;

            % copy set properties
            obj.isHRep.val = P.isHRep.val;
            obj.isVRep.val = P.isVRep.val;
            obj.emptySet.val = P.emptySet.val;
            obj.fullDim.val = P.fullDim.val;
            obj.bounded.val = P.bounded.val;
            obj.minHRep.val = P.minHRep.val;
            obj.minVRep.val = P.minVRep.val;
            return
        end

        % 2. parse input arguments: varargin -> vars
        [A,b,Ae,be,V] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(A,b,Ae,be,V,nargin);

        % 4. compute properties and hidden properties
        [A,b,Ae,be,V,isHRep,isVRep] = aux_computeProperties(A,b,Ae,be,V,nargin);
        [empty,bounded,fullDim,minHRep,minVRep,V,isHRep,isVRep] = ...
            aux_computeHiddenProperties(A,b,Ae,be,V,isHRep,isVRep);

        % 4a. assign properties
        obj.A_.val = A;
        obj.b_.val = b;
        obj.Ae_.val = Ae;
        obj.be_.val = be;
        obj.V_.val = V;

        obj.isHRep.val = isHRep;
        obj.isVRep.val = isVRep;
        obj.emptySet.val = empty;
        obj.bounded.val = bounded;
        obj.fullDim.val = fullDim;
        obj.minHRep.val = minHRep;
        obj.minVRep.val = minVRep;

        % 5. set precedence (fixed)
        obj.precedence = 80;
        
    end

    % prohibit access to the constraints and halfspaces, as they are
    % setproperty objects and thus not directly accessible as per usual

    % vertices
    function val = get.V(obj)
        if ~obj.isVRep.val
            throw(CORAerror('CORA:specialError',...
                "The vertex representation is not available. " + ...
                "Call the function 'polytope/vertices'."));
        end
        val = obj.V_.val;
    end

    % inequality constraints
    function val = get.A(obj)
        if ~obj.isHRep.val
            throw(CORAerror('CORA:specialError',...
                "The halfspace representation is not available. " + ...
                "Call the function 'polytope/constraints'."));
        end
        val = obj.A_.val;
    end
    function val = get.b(obj)
        if ~obj.isHRep.val
            throw(CORAerror('CORA:specialError',...
                "The halfspace representation is not available. " + ...
                "Call the function 'polytope/constraints'."));
        end
        val = obj.b_.val;
    end
    % equality constraints
    function val = get.Ae(obj)
        if ~obj.isHRep.val
            throw(CORAerror('CORA:specialError',...
                "The halfspace representation is not available. " + ...
                "Call the function 'polytope/constraints'."));
        end
        val = obj.Ae_.val;
    end
    function val = get.be(obj)
        if ~obj.isHRep.val
            throw(CORAerror('CORA:specialError',...
                "The halfspace representation is not available. " + ...
                "Call the function 'polytope/constraints'."));
        end
        val = obj.be_.val;
    end

end

methods (Static = true)
    P = generateRandom(varargin) % generate random polytope
    P = enclosePoints(points,varargin) % enclose point cloud with polytope
    P = empty(n) % instantiate empty polytope
    P = Inf(n) % instantiate polytope representing R^n
    P = origin(n) % instantiate polytope representing the origin in R^n
end

methods (Access = protected)
    [abbrev,printOrder] = getPrintSetInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [A,b,Ae,be,V] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

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
        if isnumeric(V)
            if any(any(isnan(V)))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Vertices have to be non-nan.'));
            elseif size(V,1) > 1 && any(any(isinf(V)))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'nD vertices for n > 1 have to be finite.'));
            end
        else
            throw(CORAerror('CORA:wrongInputInConstructor','V has to be numeric.'))
        end

        % check type
        if ~isnumeric(A)
            throw(CORAerror('CORA:wrongInputInConstructor','A has to be numeric.'))
        end
        if ~isnumeric(b)
            throw(CORAerror('CORA:wrongInputInConstructor','b has to be numeric.'))
        end
        if ~isnumeric(Ae)
            throw(CORAerror('CORA:wrongInputInConstructor','Ae has to be numeric.'))
        end
        if ~isnumeric(be)
            throw(CORAerror('CORA:wrongInputInConstructor','be has to be numeric.'))
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

function [A,b,Ae,be,V,isHRep,isVRep] = aux_computeProperties(A,b,Ae,be,V,n_in)

    % dimension
    n = max([size(A,2),size(Ae,2),size(V,1)]);

    % offsets must be column vectors
    b = reshape(b,[],1);
    be = reshape(be,[],1);

    % store which representation is given (constructor only allows one)
    isVRep = n_in == 1;
    isHRep = ~isVRep;

    % in 1D, remove redundancies (otherwise keep V as is)
    if isVRep && n == 1 && size(V,2) > 0
        V = [min(V), max(V)];
        if withinTol(V(1),V(2),eps)
            V = V(1);
        end
    end

    % empty constraint matrices must have correct dimension (necessary
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

function [empty,bounded,fullDim,minHRep,minVRep,V,isHRep,isVRep] = ...
            aux_computeHiddenProperties(A,b,Ae,be,V,isHRep,isVRep)

    % init hidden properties as unknown
    empty = []; bounded = []; fullDim = []; minHRep = []; minVRep = [];    

    % dimension
    n = max([size(A,2),size(Ae,2),size(V,1)]);

    % note: representations for 1D polytopes always minimal

    % check if instantiated via vertices
    if isVRep
        % check emptiness
        empty = size(V,2) == 0;

        % max. 1 vertex -> minimal V-representation
        minVRep = size(V,2) <= 1 || n == 1;

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
                [~,S,~] = svd(V - mean(V,2));
                fullDim = n == nnz(~withinTol(S,0,1e-12));
            end
        end

    elseif isHRep
        if isempty(A) && isempty(Ae)
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
                % only a single infeasible constraint required to represent
                % an empty set
                minHRep = length(be) + length(b) == 1;
                % init no vertices (which is the minimal representation)
                V = zeros(n,0);
                isVRep = true;
                minVRep = true;
            end
        end

    end

end

% ------------------------------ END OF CODE ------------------------------
