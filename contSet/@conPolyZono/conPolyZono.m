classdef (InferiorClasses = {?intervalMatrix, ?matZonotope})conPolyZono < contSet
% conPolyZono - class definition for constrained polynomial zonotopes 
%
% Syntax:  
%    obj = conPolyZonotope(c,G,expMat)
%    obj = conPolyZonotope(c,G,expMat,Grest)
%    obj = conPolyZonotope(c,G,expMat,Grest,id)
%    obj = conPolyZonotope(c,G,expMat,A,b,expMat_)
%    obj = conPolyZonotope(c,G,expMat,A,b,expMat_,Grest)
%    obj = conPolyZonotope(c,G,expMat,A,b,expMat_,Grest,id)
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
%
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope, zonotope, conZonotope

% Author:       Niklas Kochdumper
% Written:      06-November-2018 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c (:,1) {mustBeNumeric,mustBeFinite} = [];        % center
    G (:,:) {mustBeNumeric,mustBeFinite} = [];        % generator matrix
    expMat (:,:) {mustBeInteger} = [];                % exponent matrix
    A (:,:) {mustBeNumeric,mustBeFinite} = [];        % constraint gens.
    b (:,:) {mustBeNumeric,mustBeFinite} = [];        % constraint offset
    expMat_ (:,:) {mustBeNumeric,mustBeFinite} = [];  % constraint exponent
    Grest (:,:) {mustBeNumeric,mustBeFinite} = [];    % indep. generators
    id (:,1) {mustBeInteger} = [];                    % identifier vector
end
   
methods

    % object constructor
    function obj = conPolyZono(c,G,expMat,varargin)
        
        % parse input arguments
        Grest = []; id = []; A = []; b = []; expMat_ = [];
        
        if nargin == 1 && isa(c,'conPolyZono')
            % copy constructor
           obj = c; return;
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
        elseif nargin ~= 3
            [id,msg] = errConstructor(); error(id,msg);
        end
        
        % call superclass method
        if isempty(id)
           id = (1:size(expMat,1))'; 
        end
        
        % check correctness of user input 
        if ~all(all(floor(expMat) == expMat)) || ~all(all(expMat >= 0)) || ...
            size(expMat,2) ~= size(G,2)
            [id,msg] = errConstructor('Invalid exponent matrix.'); error(id,msg);
        end
        
        if ~isempty(expMat_)
            if size(expMat,1) ~= size(expMat_,1) 
                [id,msg] = errConstructor('Input arguments "expMat" and "expMat_" are not compatible.');
                error(id,msg);
            end
            if ~all(all(floor(expMat_) == expMat_)) || ~all(all(expMat_ >= 0))
                [id,msg] = errConstructor('Invalid constraint exponent matrix.');
                error(id,msg);
            end
            if isempty(A) || size(A,2) ~= size(expMat_,2)
                [id,msg] = errConstructor('Input arguments "A" and "expMat_" are not compatible.');
                error(id,msg);
            end
            if isempty(b) || size(b,2) > 1 || size(b,1) ~= size(A,1)
                [id,msg] = errConstructor('Input arguments "A" and "b" are not compatible.');
                error(id,msg);
            end
        elseif ~isempty(A) || ~isempty(b)
            [id,msg] = errConstructor('Invalid constraint exponent matrix.');
            error(id,msg);
        end
        
        % assign object properties
        obj.c = c; obj.G = G; obj.expMat = expMat; obj.A = A; obj.b = b;
        obj.expMat_ = expMat_; obj.Grest = Grest; obj.id = id;
        
        % set parent object properties
        obj.dimension = length(obj.c);
    end
end

methods (Static = true)
    cPZ = generateRandom(varargin) % generate random conPolyZono object
end
end
%------------- END OF CODE --------------