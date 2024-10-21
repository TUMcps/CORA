classdef zonotope < contSet
% zonotope - object constructor for zonotope objects
%
% Description:
%    This class represents zonotopes objects defined as
%    {c + \sum_{i=1}^p beta_i * g^(i) | beta_i \in [-1,1]}.
%
% Syntax:
%    obj = zonotope(c,G)
%    obj = zonotope(Z)
%
% Inputs:
%    c - center vector
%    G - generator matrix
%    Z - center vector and generator matrix Z = [c,G]
%
% Outputs:
%    obj - generated zonotope object
%
% Example: 
%    c = [1;1];
%    G = [1 1 1; 1 -1 0];
%    Z = zonotope(c,G);
%    plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope, conZonotope, zonoBundle, conPolyZono

% Authors:       Matthias Althoff, Niklas Kochdumper, Tobias Ladner
% Written:       14-September-2006 
% Last update:   22-March-2007
%                04-June-2010
%                08-February-2011
%                18-November-2015
%                05-December-2017 (DG, redefined class with the new standard)
%                28-April-2019 (code shortened)
%                01-May-2020 (NK, new constructor + removed orientation prop.)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                29-March-2023 (TL, optimized constructor)
%                13-September-2023 (TL, replaced Z property with c and G)
%                05-October-2024 (MW, removed halfspace property)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = {?contSet, ?matrixSet}, GetAccess = public)
    c;                  % center
    G;                  % generator matrix

    % legacy
    Z = [];             % concatenated center and generator matrix
    halfspace = [];     % halfspace representation
end

methods

    function obj = zonotope(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end
        assertNarginConstructor(1:2,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'zonotope')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [c,G] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(c,G,nargin);

        % 4. compute properties
        [c,G] = aux_computeProperties(c,G);

        % 5. assign properties
        obj.c = c;
        obj.G = G;

        % 6. set precedence (fixed)
        obj.precedence = 110;
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random zonotope
    Z = enclosePoints(points,varargin) % enclose point cloud with zonotope
    Z = empty(n) % instantiates an empty zonotope
    Z = origin(n);
end

methods (Access = protected)
    [abbrev,printOrder] = getPrintSetInfo(S)
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
        CORAwarning('CORA:deprecated', 'property', 'zonotope.Z', 'CORA v2024', ...
            'Please use zonotope.c and zonotope.G instead.', ...
            'This change was made to be consistent with the notation in papers.')
        Z = [obj.c, obj.G];
    end

    function obj = set.Z(obj, Z)
        CORAwarning('CORA:deprecated', 'property', 'zonotope.Z', 'CORA v2024', ...
            'Please use zonotope.c and zonotope.G instead.', ...
            'This change was made to be consistent with the notation in papers.')
        if ~isempty(Z)
            obj.c = Z(:,1);
            obj.G = Z(:,2:end);
        end
    end

    % getter for removed halfspace property
    function hs = get.halfspace(obj)
        CORAwarning('CORA:deprecated', 'property', 'zonotope.halfspace', 'CORA v2025', ...
            'Please call polytope(obj) instead.', ...
            'This change was made to avoid code redundancy.')
        hs = obj.halfspace;
    end
    function obj = set.halfspace(obj, hs)
        CORAwarning('CORA:deprecated', 'property', 'zonotope.halfspace', 'CORA v2025', ...
            'Please use polytope objects instead.', ...
            'This change was made to avoid code redundancy.')
        obj.halfspace = [];
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % no input arguments
    if nargin == 0
        c = zeros(0,0); G = zeros(0,0);
        return
    end

    % set default values depending on nargin
    if nargin == 1
        if size(varargin{1},2) == 0
            c = varargin{1};
            G = [];
        else
            c = varargin{1}(:,1);
            G = varargin{1}(:,2:end);
        end
    elseif nargin == 2
        [c,G] = setDefaultValues({[],[]},varargin);
    end

end

function aux_checkInputArgs(c,G,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        if n_in == 1

            inputArgsCheck({{[c,G], 'att', 'numeric', 'nonnan'}})

        elseif n_in == 2
        
            inputArgsCheck({ ...
                {c, 'att', 'numeric', 'nonnan'}; ...
                {G, 'att', 'numeric', 'nonnan'}; ...
            })
    
            % check dimensions
            if isempty(c) && ~isempty(G)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Center is empty.'));
            elseif ~isempty(c) && ~isvector(c)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Center is not a vector.'));
            elseif ~isempty(G) && size(c,1) ~= size(G,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Dimension mismatch between center and generator matrix.'));  
            end

        end
        
    end

end

function [c,G] = aux_computeProperties(c,G)

    % if G is empty, set correct dimension
    if isempty(G)
        G = zeros(size(c,1),0);
    end

end

% ------------------------------ END OF CODE ------------------------------
