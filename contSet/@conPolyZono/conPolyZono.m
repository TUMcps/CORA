classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) conPolyZono < contSet
% conPolyZono - class definition for constrained polynomial zonotopes 
%
% Syntax:
%    obj = conPolyZono(c)
%    obj = conPolyZono(c,G,E)
%    obj = conPolyZono(c,G,E,GI)
%    obj = conPolyZono(c,G,E,GI,id)
%    obj = conPolyZono(c,G,E,A,b,EC)
%    obj = conPolyZono(c,G,E,A,b,EC,GI)
%    obj = conPolyZono(c,G,E,A,b,EC,GI,id)
%
% Inputs:
%    c - constant offset (dimension: [n,1])
%    G - generator matrix (dimension: [n,h])
%    E - exponent matrix (dimension: [p,h])
%    A - constraint generator matrix (dimension: [m,q])
%    b - constraint offset (dimension: [m,1])
%    EC - constraint exponent matrix (dimension: [p,q])
%    GI - generator matrix of independent generators (dimension: [n,v])
%    id - vector containing the integer identifiers for all factors
%         (dimension: [p,1])
%
% Outputs:
%    obj - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 1 1 1];
%    E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [0 1 2; 1 0 0; 0 1 0];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%    plot(cPZ,[1,2],'FaceColor','r','Splits',8);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope, zonotope, conZonotope

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       06-November-2018 
% Last update:   14-December-2022 (TL, property check in inputArgsCheck)
%                13-September-2023 (TL, replaced Grest/expMat/_ properties with GI/E/EC)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    c;        % center
    G;        % generator matrix
    E;        % exponent matrix
    A;        % constraint generators
    b;        % constraint offset
    EC;       % constraint exponent matrix
    GI;       % independent generators
    id;       % identifier vector

    % legacy
    expMat;
    expMat_;
    Grest;
end
   
methods
    % object constructor
    function obj = conPolyZono(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'conPolyZono')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [c,G,E,A,b,EC,GI,id] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(c,G,E,A,b,EC,GI,id,nargin);

        % 4. compute properties
        [c,G,E,A,b,EC,GI,id] = aux_computeProperties(c,G,E,A,b,EC,GI,id);

        % 5. assign properties
        obj.c = c; obj.G = G; obj.E = E; obj.A = A; obj.b = b;
        obj.EC = EC; obj.GI = GI; obj.id = id;
    end
end

methods (Static = true)
    cPZ = generateRandom(varargin) % generate random conPolyZono object
    C = empty(n) % instantiates an empty constrained polynomial zonotope
end


% getter & setter ---------------------------------------------------------

methods
    % getter & setter for legacy Grest property
    function Grest = get.Grest(obj)
        warning(['CORA: The property conPolyZono.Grest is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conPolyZono.GI instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        Grest = obj.GI;
    end

    function obj = set.Grest(obj, Grest)
        warning(['CORA: The property conPolyZono.Grest is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conPolyZono.GI instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        obj.GI = Grest;
    end

    % getter & setter for legacy expMat property
    function expMat = get.expMat(obj)
        warning(['CORA: The property conPolyZono.expMat is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conPolyZono.E instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        expMat = obj.E;
    end

    function obj = set.expMat(obj, expMat)
        warning(['CORA: The property conPolyZono.expMat is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conPolyZono.E instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        obj.E = expMat;
    end

    % getter & setter for legacy expMat_ property
    function expMat_ = get.expMat_(obj)
        warning(['CORA: The property conPolyZono.expMat_ is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conPolyZono.EC instead.']);

        expMat_ = obj.EC;
    end

    function obj = set.expMat_(obj, expMat_)
        warning(['CORA: The property conPolyZono.expMat_ is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use conPolyZono.EC instead.']);

        obj.EC = expMat_;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G,E,A,b,EC,GI,id] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 8
        % too many input arguments
        throw(CORAerror('CORA:tooManyInputArgs',8));
    end

    % default values
    [c,G,E] = setDefaultValues({[],[],[]},varargin);
    A = []; b = []; EC = []; GI = []; id = [];
    
    % set other (optional) values
    if nargin == 1 && isnumeric(c) && isvector(c)
        return;
    elseif nargin == 4 || nargin == 5
        [GI,id] = setDefaultValues({GI,id},varargin(4:end));
    elseif nargin >= 6
        [A,b,EC,GI,id] = setDefaultValues({A,b,EC,GI,id},varargin(4:end));
    end
    
    % set identifiers
    if isempty(id)
        id = (1:size(E,1))'; 
    end

end

function aux_checkInputArgs(c,G,E,A,b,EC,GI,id,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % check correctness of user input 
        inputChecks = { ...
            {c, 'att', 'numeric', {'finite'}}; ...
            {G, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {E, 'att', 'numeric', {'integer', 'matrix'}}; ...
        };
        if n_in > 5
            % only add constraints checks if they were in the input
            % to correctly indicate the position of the wrong input
            inputChecks = [
                inputChecks;
                {{A, 'att', 'numeric', {'finite', 'matrix'}}; ...
                {b, 'att', 'numeric', {'finite', 'matrix'}}; ...
                {EC, 'att', 'numeric', {'finite', 'matrix'}}}; ...
            ];
        end
        inputChecks = [
            inputChecks
            {{GI, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {id, 'att', 'numeric', {'finite'}}};
        ];
        inputArgsCheck(inputChecks)

        % center must be a vector
        if isempty(c)
            if ~isempty(G) || ~isempty(E) || ~isempty(A) || ~isempty(b) ...
                    || ~isempty(EC) || ~isempty(GI) || ~isempty(id)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Either all or none input arguments are empty.'));
            end
        elseif ~isvector(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Center must be a vector.'));
        end

        % check dimensions
        if size(E,2) ~= size(G,2)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Invalid exponent matrix.'));
        end
        if ~isempty(EC)
            if size(E,1) ~= size(EC,1) 
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Input arguments "E" and "EC" are not compatible.'));
            end
            if ~all(all(floor(EC) == EC)) || ~all(all(EC >= 0))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Invalid constraint exponent matrix.'));
            end
            if isempty(A) || size(A,2) ~= size(EC,2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Input arguments "A" and "EC" are not compatible.'));
            end
            if isempty(b) || size(b,2) > 1 || size(b,1) ~= size(A,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Input arguments "A" and "b" are not compatible.'));
            end
        elseif ~isempty(A) || ~isempty(b)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Invalid constraint exponent matrix.'));
        end
        
    end

end

function [c,G,E,A,b,EC,GI,id] = aux_computeProperties(c,G,E,A,b,EC,GI,id)

    if ~isempty(c)
        % make center a column vector
        c = reshape(c,[],1);
    else
        % set generator matrices to correct dimensions
        n = size(c,1);
        if isempty(G)
            G = zeros(n,0);
        end
        if isempty(GI)
            GI = zeros(n,0);
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
