classdef (InferiorClasses = {?intervalMatrix, ?matZonotope})polyZonotope
% polyZonotope - Object Constructor for polonomial zonotopes 
%
% Syntax:  
%    object constructor:    Obj = polyZonotope(c,G,Grest)
%                           Obj = polyZonotope(c,G,Grest,expMat)
%                           Obj = polyZonotope(c,G,Grest,expMat,id)
%
% Inputs:
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
%    Obj - Polynomial Zonotope Object
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

% Author:        Niklas Kochdumper
% Written:       26-March-2018 
% Last update:   02-May-2020 (MW, add property validation, def constructor)
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

    function Obj = polyZonotope(c,G,Grest,varargin)
        
        if nargin == 1 
            % Copy Constructor
            if isa(c,'polyZonotope')
                Obj = c;
                
            % Single point
            elseif isvector(c)
                Obj.c = c;
            end
        
        % No exponent matrix provided
        elseif nargin == 2 || nargin == 3
            
            Obj.c = c;
            Obj.G = G;
            if nargin == 3
                Obj.Grest = Grest;
            end
            
            % construct exponent matrix under the assumption
            % that all generators are independent
            Obj.expMat = eye(size(G,2));
            Obj.id = (1:size(G,2))';
            
        
        % Exponent matrix as user input
        elseif nargin == 4 || nargin == 5
            
            expMat = varargin{1};
            
            % check correctness of user input
            if ~all(all(floor(expMat) == expMat)) || ~all(all(expMat >= 0)) || ...
               size(expMat,2) ~= size(G,2)
                error('Invalid exponent matrix!');
            end
            
            % remove redundant exponents
            if ~isempty(expMat)
                [expMat,G] = removeRedundantExponents(expMat,G);
            end
            
            % assign properties
            Obj.c = c;
            Obj.G = G;
            Obj.Grest = Grest;
            Obj.expMat = expMat;
            
            % vector of integer identifiers
            if nargin == 5
                
                id = varargin{2};
                
                % check for correctness
                if length(id) ~= size(expMat,1)
                   error('Invalid vector of identifiers!'); 
                end
                
                Obj.id = id;
                
            else
                Obj.id = (1:size(expMat,1))';
            end
        

        % Else if not enough or too many inputs are passed    
        elseif nargin > 5
            disp('This class needs less input values.');
            Obj=[];
        end
        
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random polyZonotope
end


end
%------------- END OF CODE --------------