classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) conPolyZono < contSet
% conPolyZono - class definition for constrained polynomial zonotopes 
%
% Syntax:  
%    obj = conPolyZono(c)
%    obj = conPolyZono(c,G,expMat)
%    obj = conPolyZono(c,G,expMat,Grest)
%    obj = conPolyZono(c,G,expMat,Grest,id)
%    obj = conPolyZono(c,G,expMat,A,b,expMat_)
%    obj = conPolyZono(c,G,expMat,A,b,expMat_,Grest)
%    obj = conPolyZono(c,G,expMat,A,b,expMat_,Grest,id)
%
% Inputs:
%    c - constant offset (dimension: [n,1])
%    G - generator matrix (dimension: [n,h])
%    expMat - exponent matrix (dimension: [p,h])
%    A - constraint generator matrix (dimension: [m,q])
%    b - constraint offset (dimension: [m,1])
%    expMat_ - constraint exponent matrix (dimension: [p,q])
%    Grest - generator matrix of independent generators (dimension: [n,v])
%    id - vector containing the integer identifiers for all factors
%         (dimension: [p,1])
%
% Outputs:
%    obj - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 1 1 1];
%    expMat = [1 0 1 2; 0 1 1 0; 0 0 1 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    expMat_ = [0 1 2; 1 0 0; 0 1 0];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%    plot(cPZ,[1,2],'FaceColor','r','Splits',8);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope, zonotope, conZonotope

% Author:       Niklas Kochdumper
% Written:      06-November-2018 
% Last update:  14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:16-June-2023 (MW, restructure using auxiliary functions)

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c;        % center
    G;        % generator matrix
    expMat;   % exponent matrix
    A;        % constraint generators
    b;        % constraint offset
    expMat_;  % constraint exponent
    Grest;    % independent generators
    id;       % identifier vector
end
   
methods
    % object constructor
    function obj = conPolyZono(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'conPolyZono')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [c,G,expMat,A,b,expMat_,Grest,id] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(c,G,expMat,A,b,expMat_,Grest,id,nargin);

        % 4. assign properties
        obj.c = c; obj.G = G; obj.expMat = expMat; obj.A = A; obj.b = b;
        obj.expMat_ = expMat_; obj.Grest = Grest; obj.id = id;
    end
end

methods (Static = true)
    cPZ = generateRandom(varargin) % generate random conPolyZono object
end

end


% Auxiliary Functions -----------------------------------------------------

function [c,G,expMat,A,b,expMat_,Grest,id] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 8
        % too many input arguments
        throw(CORAerror('CORA:tooManyInputArgs',8));
    end

    % default values
    [c,G,expMat] = setDefaultValues({[],[],[]},varargin);
    A = []; b = []; expMat_ = []; Grest = []; id = [];
    
    % set other (optional) values
    if nargin == 0 || (nargin == 1 && isnumeric(c) && isvector(c))
        return;
    elseif nargin == 4 || nargin == 5
        [Grest,id] = setDefaultValues({Grest,id},varargin(4:end));
    elseif nargin >= 6
        [A,b,expMat_,Grest,id] = setDefaultValues({A,b,expMat_,Grest,id},varargin(4:end));
    end
    
    % set identifiers
    if isempty(id)
        id = (1:size(expMat,1))'; 
    end

end

function aux_checkInputArgs(c,G,expMat,A,b,expMat_,Grest,id,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % check correctness of user input 
        inputChecks = { ...
            {c, 'att', 'numeric', {'finite', 'column'}}; ...
            {G, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {expMat, 'att', 'numeric', {'integer', 'matrix'}}; ...
        };
        if n_in > 5
            % only add constraints checks if they were in the input
            % to correctly indicate the position of the wrong input
            inputChecks = [
                inputChecks;
                {{A, 'att', 'numeric', {'finite', 'matrix'}}; ...
                {b, 'att', 'numeric', {'finite', 'matrix'}}; ...
                {expMat_, 'att', 'numeric', {'finite', 'matrix'}}}; ...
            ];
        end
        inputChecks = [
            inputChecks
            {{Grest, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {id, 'att', 'numeric', {'finite'}}};
        ];
        inputArgsCheck(inputChecks)


        % check dimensions
        if size(expMat,2) ~= size(G,2)
             throw(CORAerror('CORA:wrongInputInConstructor',...
                 'Invalid exponent matrix.'));
        end
        if ~isempty(expMat_)
            if size(expMat,1) ~= size(expMat_,1) 
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Input arguments "expMat" and "expMat_" are not compatible.'));
            end
            if ~all(all(floor(expMat_) == expMat_)) || ~all(all(expMat_ >= 0))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Invalid constraint exponent matrix.'));
            end
            if isempty(A) || size(A,2) ~= size(expMat_,2)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Input arguments "A" and "expMat_" are not compatible.'));
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

%------------- END OF CODE --------------