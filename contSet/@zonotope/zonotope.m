classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonotope < contSet
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
% See also: ---

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
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    c, G; % zonotope center and generator 
    Z; % legacy Z = [c,g_1,...,g_p]

    % internally-set properties
    halfspace;  % halfspace representation of the zonotope
end

methods

    function obj = zonotope(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

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
        obj.halfspace = [];
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random zonotope
    Z = enclosePoints(points,varargin) % enclose point cloud with zonotope
    Z = empty(n) % instantiates an empty zonotope
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
        warning(['CORA: The property zonotope.Z is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use zonotope.c and zonotope.G instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        Z = [obj.c, obj.G];
    end

    function obj = set.Z(obj, Z)
        warning(['CORA: The property zonotope.Z is deprecated (since CORA 2024) and will be removed in a future release. ' ...
            'Please use zonotope.c and zonotope.G instead. ' ...
            'This change was made to be consistent with the notation in papers.']);

        if ~isempty(Z)
            obj.c = Z(:,1);
            obj.G = Z(:,2:end);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 2
        throw(CORAerror('CORA:tooManyInputArgs',2));
    end

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
