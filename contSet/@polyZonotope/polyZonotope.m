classdef polyZonotope < contSet
% polyZonotope - object constructor for polynomial zonotopes
%
% Definition: see CORA manual, Sec. 2.2.1.5.
%
% Syntax:
%    obj = polyZonotope(pZ)
%    obj = polyZonotope(c)
%    obj = polyZonotope(c,G)
%    obj = polyZonotope(c,G,GI)
%    obj = polyZonotope(c,[],GI)
%    obj = polyZonotope(c,G,GI,E)
%    obj = polyZonotope(c,G,[],E)
%    obj = polyZonotope(c,G,GI,E,id)
%    obj = polyZonotope(c,G,[],E,id)
%
% Inputs:
%    pZ - polyZonotope object
%    c - center of the polynomial zonotope (dimension: [nx,1])
%    G - generator matrix containing the dependent generators 
%       (dimension: [nx,N])
%    GI - generator matrix containing the independent generators
%            (dimension: [nx,M])
%    E - matrix containing the exponents for the dependent generators
%             (dimension: [p,N])
%    id - vector containing the integer identifiers for the dependent
%         factors (dimension: [p,1])
%
% Outputs:
%    obj - polyZonotope object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    GI = [0;0.5];
%    E = [1 0 3;0 1 1];
% 
%    pZ = polyZonotope(c,G,GI,E);
% 
%    plot(pZ,[1,2],'FaceColor','r');
%
% References:
%    [1] Kochdumper, N., et al. (2020). Sparse polynomial zonotopes: A novel
%        set representation for reachability analysis. IEEE Transactions on 
%        Automatic Control.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Tobias Ladner
% Written:       26-March-2018 
% Last update:   02-May-2020 (MW, add property validation, def constructor)
%                21-March-2021 (MW, error messages, size checks, restructuring)
%                14-December-2022 (TL, restructuring)
%                29-March-2023 (TL, optimized constructor)
%                13-September-2023 (TL, replaced Grest/expMat properties with GI/E)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = {?contSet, ?matrixSet}, GetAccess = public)
    c;      % center
    G;      % dependent generator matrix
    GI;     % independent generator matrix
    E;      % exponent matrix
    id;     % identifier vector

    % legacy
    Grest = [];
    expMat = [];
end
   
methods

    function obj = polyZonotope(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end
        assertNarginConstructor(1:5,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'polyZonotope')
            obj = varargin{1}; return
        end

        % 1. parse input arguments: varargin -> vars
        [c,G,GI,E,id] = aux_parseInputArgs(varargin{:});

        % 2. check correctness of input arguments
        aux_checkInputArgs(c,G,GI,E,id,nargin);

        % 3. compute properties
        [c,G,GI,E,id] = aux_computeProperties(c,G,GI,E,id);

        % 4. assign properties
        obj.c = c;
        obj.G = G;
        obj.GI = GI;
        obj.E = E;
        obj.id = id;

        % 5. set precedence (fixed)
        obj.precedence = 70;

    end
end

methods (Static = true)
    pZ = generateRandom(varargin) % generate random polyZonotope
    pZ = empty(n) % instantiates an empty polyZonotope
    pZ = origin(n) % instantiates a polyZonotope representing the origin
end

methods (Access = protected)
    [abbrev,printOrder] = getPrintSetInfo(S)
end

% plotting
methods (Access = {?contSet})
    han = plot3D(S,varargin);
end


% getter & setter ---------------------------------------------------------

methods
    % getter & setter for legacy Grest property
    function Grest = get.Grest(obj)
        CORAwarning('CORA:deprecated', 'property', 'polyZonotope.Grest', 'CORA v2024', ...
            'Please use polyZonotope.GI instead.', ...
            'This change was made to be consistent with the notation in papers.')
        Grest = obj.GI;
    end

    function obj = set.Grest(obj, Grest)
        CORAwarning('CORA:deprecated', 'property', 'polyZonotope.Grest', 'CORA v2024', ...
            'Please use polyZonotope.GI instead.', ...
            'This change was made to be consistent with the notation in papers.')
        obj.GI = Grest;
    end

    % getter & setter for legacy expMat property
    function expMat = get.expMat(obj)
        CORAwarning('CORA:deprecated', 'property', 'polyZonotope.expMat', 'CORA v2024', ...
            'Please use polyZonotope.E instead.', ...
            'This change was made to be consistent with the notation in papers.')
        expMat = obj.E;
    end

    function obj = set.expMat(obj, expMat)
        CORAwarning('CORA:deprecated', 'property', 'polyZonotope.expMat', 'CORA v2024', ...
            'Please use polyZonotope.E instead.', ...
            'This change was made to be consistent with the notation in papers.')
        obj.E = expMat;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G,GI,E,id] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % no input arguments
    if nargin == 0
        c = []; G = []; GI = []; E = []; id = []; 
        return
    end

    % set default values
    [c,G,GI,E,id] = setDefaultValues({[],[],[],[],[]},varargin);

end

function aux_checkInputArgs(c,G,GI,E,id,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({ ...
            {c, 'att', 'numeric', 'finite'};
            {G, 'att', 'numeric', {'finite', 'matrix'}};
            {GI, 'att', 'numeric', {'finite', 'matrix'}};
            {E, 'att', 'numeric', ...
                {'integer', 'nonnegative', 'matrix'}};
            {id, 'att', 'numeric', 'integer'};
        })

        % check dimensions ---
        % c
        if ~isempty(c) && size(c, 2) ~= 1
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Center should be a column vector.'));
        end
        % G
        if ~isempty(G) && size(G,1) ~= size(c,1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Dimension mismatch between center and dependent generator matrix.'));
        end
        % GI
        if ~isempty(GI) && size(GI,1) ~= size(c,1)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Dimension mismatch between center and independent generator matrix.'));
        end
        % E
        if ~isempty(E) && size(E,2) ~= size(G,2)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Dimension mismatch between dependent generator matrix and exponent matrix.'));
        end
        % id
        if ~isempty(id) && size(id, 2) ~= 1
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Identifier vector should be a column vector.'));
        end
        if ~isempty(E) && ~isempty(id) && size(id, 1) ~= size(E,1)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Dimension mismatch between exponent matrix and identifier vector.'));
        end
    end

end

function [c,G,GI,E,id] = aux_computeProperties(c,G,GI,E,id)

    % remove redundancies
    if ~isempty(E)
        [E,G] = removeRedundantExponents(E,G);
    end

    % if G/GI is empty, set correct dimension
    if isempty(G)
        G = zeros(size(c,1),0);
    end
    if isempty(GI)
        GI = zeros(size(c,1),0);
    end

    % default value for exponent matrix
    if isempty(E)
        E = eye(size(G,2));
    end

    % number of dependent factors
    if isempty(id)
        p = size(E,1);
        id = (1:p)'; % column vector
    end

end

% ------------------------------ END OF CODE ------------------------------
