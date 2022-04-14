classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) polyZonotope < contSet
% polyZonotope - object constructor for polynomial zonotope object
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
%    obj - polynomial zonotope Object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    Grest = [0;0.5];
%    expMat = [1 0 3;0 1 1];
%
%    pZ = polyZonotope(c,G,Grest,expMat);
%
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');  
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:        Niklas Kochdumper, Mark Wetzlinger
% Written:       26-March-2018 
% Last update:   02-May-2020 (MW, add property validation, def constructor)
%                21-March-2021 (MW, errConstructor, size checks, restructuring)
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c (:,1) {mustBeNumeric,mustBeFinite} = [];
    G (:,:) {mustBeNumeric,mustBeFinite} = [];
    Grest (:,:) {mustBeNumeric,mustBeFinite} = [];
    expMat (:,:) {mustBeInteger} = [];
    id (:,1) {mustBeInteger} = [];
end
   
methods

    function obj = polyZonotope(c,G,Grest,expMat,id)
        
        if nargin == 0
            % return empty set
            obj.c = [];
            c = []; % for contSet / dimension
        
        elseif nargin == 1 
            % Copy Constructor
            if isa(c,'polyZonotope')
                obj = c;
            % Single point
            else
                obj.c = c;
            end
        
        elseif nargin >= 2 && nargin <= 5
            
            % check center for emptyness
            if isempty(c)
                [id,msg] = errConstructor('Center is empty.');
                error(id,msg);
            end
            % assign center
            obj.c = c;
            
            % check sizes of c and G
            if ~isempty(G) && length(c) ~= size(G,1)
                [id,msg] = errConstructor(...
                    'Dimension mismatch between center and dependent generator matrix.');
                error(id,msg);
            end
            % assign dependent generator matrix
            % note: later overwritten if expMat given (to avoid checking twice)
            obj.G = G;
            
            if nargin >= 3
                % check sizes of c and Grest
                if ~isempty(Grest) && length(c) ~= size(Grest,1)
                    [id,msg] = errConstructor(...
                        'Dimension mismatch between center and independent generator matrix.');
                    error(id,msg);
                end
                % assign independent generator matrix
                obj.Grest = Grest;
            end
            
            if nargin <= 3
                % construct exponent matrix and identifiers under
                % the assumption that all generators are independent
                obj.expMat = eye(size(G,2));
                obj.id = (1:size(G,2))';
            
            else

                % check correctness of user input
                if ~all(all(floor(expMat) == expMat)) || ~all(all(expMat >= 0)) ...
                        || size(expMat,2) ~= size(G,2)
                    [id,msg] = errConstructor('Invalid exponent matrix.');
                    error(id,msg);
                end

                % remove redundant exponents
                if ~isempty(expMat)
                    [expMat,G] = removeRedundantExponents(expMat,G);
                end

                % assign properties
                obj.G = G; % this overwrites the former assignment
                obj.expMat = expMat;

                % vector of integer identifiers
                if nargin == 5

                    % check for correctness
                    if length(id) ~= size(expMat,1)
                        [errid,msg] = errConstructor('Invalid vector of identifiers.');
                        error(errid,msg);
                    end
                    % assign value
                    obj.id = id;

                else
                    % no identifiers provided
                    obj.id = (1:size(expMat,1))';
                end
                
            end

        % Too many inputs are passed    
        elseif nargin > 5
            
            % too many input arguments
            [id,msg] = errConstructor('Too many input arguments.');
            error(id,msg);
        end
        
        % set parent object properties
        obj.dimension = length(c);
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random polyZonotope
end


end
%------------- END OF CODE --------------