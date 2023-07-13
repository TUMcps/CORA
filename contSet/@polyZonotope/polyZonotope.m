classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) polyZonotope < contSet
% polyZonotope - object constructor for polynomial zonotopes
%
% Definition: see CORA manual, Sec. 2.2.1.5.
%
% Syntax:
%    obj = polyZonotope()
%    obj = polyZonotope(pZ)
%    obj = polyZonotope(c,G)
%    obj = polyZonotope(c,G,Grest)
%    obj = polyZonotope(c,[],Grest)
%    obj = polyZonotope(c,G,Grest,expMat)
%    obj = polyZonotope(c,G,[],expMat)
%    obj = polyZonotope(c,G,Grest,expMat,id)
%    obj = polyZonotope(c,G,[],expMat,id)
%
% Inputs:
%    pZ - polyZonotope object
%    c - center of the polynomial zonotope (dimension: [nx,1])
%    G - generator matrix containing the dependent generators 
%       (dimension: [nx,N])
%    Grest - generator matrix containing the independent generators
%            (dimension: [nx,M])
%    expMat - matrix containing the exponents for the dependent generators
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
%    Grest = [0;0.5];
%    expMat = [1 0 3;0 1 1];
% 
%    pZ = polyZonotope(c,G,Grest,expMat);
% 
%    plot(pZ,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:        Niklas Kochdumper, Mark Wetzlinger, Tobias Ladner
% Written:       26-March-2018 
% Last update:   02-May-2020 (MW, add property validation, def constructor)
%                21-March-2021 (MW, error messages, size checks, restructuring)
%                14-December-2022 (TL, restructuring)
%                29-March-2023 (TL: optimized constructor)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c = [];
    G = [];
    Grest = [];
    expMat = [];
    id = [];
end
   
methods

    function obj = polyZonotope(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'polyZonotope')
            obj = varargin{1}; return
        end

        % 1. parse input arguments: varargin -> vars
        [c,G,Grest,expMat,id] = aux_parseInputArgs(varargin{:});

        % 2. check correctness of input arguments
        aux_checkInputArgs(c,G,Grest,expMat,id,nargin);

        % 3. remove redundancies
        if ~isempty(expMat)
            [expMat,G] = removeRedundantExponents(expMat,G);
        end

        % 4. assign properties
        obj.c = c;
        obj.G = G;
        obj.Grest = Grest;
        obj.expMat = expMat;
        obj.id = id;

    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random polyZonotope
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G,Grest,expMat,id] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 5
        throw(CORAerror('CORA:tooManyInputArgs',5));
    end

    % no input arguments
    if nargin == 0
        c = []; G = []; Grest = []; expMat = []; id = []; 
        return
    end

    % set default values
    [c,G,Grest,expMat,id] = setDefaultValues({[],[],[],[],[]},varargin);

    % default value for exponent matrix
    if all(size(expMat) == [0,0])
        h = size(G,2);
        expMat = eye(h);
    end

    % number of dependent factors
    if isempty(id)
        p = size(expMat,1);
        id = (1:p)'; % column vector
    end

end

function aux_checkInputArgs(c,G,Grest,expMat,id,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({ ...
            {c, 'att', 'double', 'finite'};
            {G, 'att', 'double', {'finite', 'matrix'}};
            {Grest, 'att', 'double', {'finite', 'matrix'}};
            {expMat, 'att', 'double', ...
                {'integer', 'nonnegative', 'matrix'}};
            {id, 'att', 'double', 'integer'};
        })

        % check dimensions
        if ~isempty(c) && size(c, 2) ~= 1
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Center should be a column vector.'));
        end
        if size(id, 2) ~= 1
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Identifier vector should be a column vector.'));
        end
        if ~isempty(G) && size(G,1) ~= size(c,1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Dimension mismatch between center and dependent generator matrix.'));
        end
        if ~isempty(Grest) && size(Grest,1) ~= size(c,1)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Dimension mismatch between center and independent generator matrix.'));
        end
        if size(expMat,2) ~= size(G,2)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Dimension mismatch between dependent generator matrix and exponent matrix.'));
        end
        if size(id, 1) ~= size(expMat,1)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Dimension mismatch between exponent matrix and identifier vector.'));
        end
    end

end

%------------- END OF CODE --------------