classdef (InferiorClasses = {?mp}) matZonotope
% matZonotope class 
%
% Syntax:
%    obj = matZonotope()
%    obj = matZonotope(C,G)
%
% Inputs:
%    C - center matrix (n x m)
%    G - h generator matrices stored as (n x m x h)
%
% Outputs:
%    obj - generated matZonotope object
%
% Example:
%    C = [0 0; 0 0];
%    G(:,:,1) = [1 3; -1 2];
%    G(:,:,2) = [2 0; 1 -1];
%
%    matZ = matZonotope(C,G);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix, matPolytope

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       14-September-2006 
% Last update:   22-March-2007
%                04-June-2010
%                27-August-2019
%                03-April-2023 (MW, remove property dim)
%                25-April-2024 (TL, matZ.C & matZ.G property, speed up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    C               % center (n x m)
    G               % generators (n x m x h)

    % legacy
    gens
    generator = []; % was stored as cell array before
end
    
methods
    
    % class constructor
    function obj = matZonotope(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'matZonotope')
            obj = varargin{1}; return
        end

        % 2. parse input arguments
        [C,G] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(C,G)

        % 4. update properties
        [C,G] = aux_computeProperties(C,G);

        % 5. assign properties
        obj.C = C;
        obj.G = G;
    end
         
    %methods in seperate files     
    matZ = plus(summand1,summand2)
    matZ = mtimes(factor1,factor2)
    matZ = mpower(matZ,exponent)
    matZ = powers(matZ,varargin)
    matZ = expmInd(matZ,maxOrder)
    [eZ,eI,zPow,iPow,E] = expmMixed(matZ,r,intermediateOrder,maxOrder)
    [eZ,eI,zPow,iPow,E] = expmIndMixed(matZ,intermediateOrder,maxOrder)
    [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,u)
    intMat = intervalMatrix(matZ,varargin)
    matZ = zonotope(matZ)
    dist = expmDist(matZ,intMat,maxOrder)
    matZred = reduce(matZ,option,order,filterLength)
    vol = volume(matZ)
    matZ1 = concatenate(matZ1,matZ2)
    res = norm(matZ,varargin)
    newObj = subsref(matZ,S)
    h = numgens(matZ)
    matZ_transposed = transpose(matZ)
    matZ_reshaped = reshape(matZ,n,m)
    res = representsa(matZ,varargin)
        
    %display functions
    plot(matZ,varargin)
    display(matZ)

    % legacy ---

    function Gs = get.generator(matZ)
        CORAwarning("CORA:deprecated","property",'matZonotope.generator','CORA v2024.2,0','With appropriate changes, please use matZonotope.G instead.','This change was made to improve speed.');
        Gs = matZ.G;
    end

    function matZ = set.generator(matZ,generator)
        CORAerror("CORA:noops",'matZ.generator');
    end
    
    function h = get.gens(matZ)
        CORAwarning("CORA:deprecated","property",'matZonotope.gens','CORA v2024.2,0','Please use matZonotope/numgens() instead.','This change was made to reduce internal maintenance.');
        h = numgens(matZ);
    end

    function matZ = set.gens(matZ,gens)
        CORAerror("CORA:noops",'matZ.gens');
    end

end

end


% Auxiliary functions -----------------------------------------------------

function [C,G] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

% default values
G = [];

% parse input
if nargin == 0
    C = zeros(0,0); 
elseif nargin == 1
    C = varargin{1};
elseif nargin == 2
    [C,G] = varargin{:};
elseif nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% fix generators to allow [] for no generators
if isempty(G)
    G = zeros([size(C),0]);
end
    
end

function aux_checkInputArgs(C,G)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED

        if iscell(G)
            % legacy
            inputArgsCheck({ ...
                {C, 'att', 'numeric','nonnan'}; ...
                {G, 'att', 'cell'}; ... % cell is legacy
            })

            % input checks are less rigorous here ..
    
            % check dimensions
            if isempty(C) && ~isempty(G)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Center is empty.'));
            end
            for i = 1:numel(G)
                % check each generator matrix
                Gi = G{i};
                if ~all(size(C,1:2) == size(Gi,1:2))
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        sprintf('Dimension mismatch between center and generator %i.',i)));  
                end
            end

        else
            inputArgsCheck({ ...
                {C, 'att', 'numeric','nonnan'}; ...
                {G, 'att', 'numeric','nonnan'}; ...
            })
    
            % check dimensions
            if isempty(C) && ~isempty(G)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Center is empty.'));
            end
            if ~all(size(C,1:2) == size(G,1:2))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Dimension mismatch between center and generators.'));  
            end
        end
        
    end
end

function [C,G] = aux_computeProperties(C,G)
    if iscell(G)
        % legacy, convert to (n x m x h) shape

        % show warning
        CORAwarning('CORA:deprecated','constructor for matZonotope using a','cell of generators','CORA v2024.2,0','Please use a single numeric matrix with dimensions (n x m x h) instead.','This change was made to improve speed.')

        % store given generators
        G_legacy = G;

        % preallocate new generators
        G = zeros([size(C),numel(G_legacy)]);

        % copy generators
        for i=1:numel(G_legacy)
            G(:,:,i) = G_legacy{i};
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
