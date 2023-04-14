classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) conPolyZono < contSet
% conPolyZono - class definition for constrained polynomial zonotopes 
%
% Syntax:  
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
% Last revision: ---

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
    function obj = conPolyZono(c,G,expMat,varargin)
        
        % parse input arguments
        if nargin == 2 % at least 3; 0 = empty set, 1 = copy
            % not enough input arguments
            throw(CORAerror('CORA:notEnoughInputArgs',3));
        elseif nargin > 8
            % too many input arguments
            throw(CORAerror('CORA:tooManyInputArgs',8));
        end

        if nargin == 0
            c = []; G = []; expMat = [];
        end
    
        Grest = []; id = []; A = []; b = []; expMat_ = [];
        
        if nargin == 0
            % empty set
            return
        elseif nargin == 1 && isa(c,'conPolyZono')
            % copy constructor
            obj = c; return;
        elseif nargin == 1 && isnumeric(c) && isvector(c)
            obj.c = c; return;
        elseif nargin == 4
           Grest = varargin{1};
        elseif nargin == 5
           Grest = varargin{1};
           id = varargin{2};
        elseif nargin == 6
           A = varargin{1};
           b = varargin{2};
           expMat_ = varargin{3};
        elseif nargin == 7
           A = varargin{1};
           b = varargin{2};
           expMat_ = varargin{3};
           Grest = varargin{4};
        elseif nargin == 8
           A = varargin{1};
           b = varargin{2};
           expMat_ = varargin{3};
           Grest = varargin{4};
           id = varargin{5};
        end
        
        % call superclass method
        if isempty(id)
            id = (1:size(expMat,1))'; 
        end


        % check correctness of user input 
        inputChecks = { ...
            {c, 'att', 'numeric', {'finite', 'column'}}; ...
            {G, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {expMat, 'att', 'numeric', {'integer', 'matrix'}}; ...
        };
        if nargin > 5
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
        
        % assign object properties
        obj.c = c; obj.G = G; obj.expMat = expMat; obj.A = A; obj.b = b;
        obj.expMat_ = expMat_; obj.Grest = Grest; obj.id = id;
    end
end

methods (Static = true)
    cPZ = generateRandom(varargin) % generate random conPolyZono object
end

end

%------------- END OF CODE --------------