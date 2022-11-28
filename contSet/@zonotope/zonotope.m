classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonotope < contSet
% zonotope - object constructor for zonotope objects
%
% Description:
%    This class represents zonotopes objects defined as
%    {c + \sum_{i=1}^p beta_i * g^(i) | beta_i \in [-1,1]}.
%
% Syntax:
%    obj = zonotope()
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      14-September-2006 
% Last update:  22-March-2007
%               04-June-2010
%               08-February-2011
%               18-November-2015
%               05-December-2017 (DG) class is redefined in compliance with
%               the new standard.
%               28-April-2019 code shortened
%               1-May-2020 (NK) new constructor + removed orientation prop.
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % zonotope center and generator Z = [c,g_1,...,g_p]
    Z (:,:) {mustBeNumeric,mustBeNonNan} = [];

    % halfspace representation of the zonotope
    halfspace = [];
end

methods

    function obj = zonotope(varargin)
        
        % If no argument is passed (default constructor)
        if nargin == 0
            obj.Z = [];
            obj.halfspace = [];

        % If 1 argument is passed
        else
            
            if nargin == 1
            
                % input is a zonotope -> copy object
                if isa(varargin{1},'zonotope')
                    obj = varargin{1};
                else
                    % List elements of the class
                    obj.Z = varargin{1}; 
                    obj.halfspace = [];
                end

            % If 2 arguments are passed
            elseif nargin == 2

                % wrong inputs
                if isempty(varargin{1})
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Center is empty.'));
                elseif ~isvector(varargin{1})
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Center is not a vector.'));
                elseif ~isempty(varargin{2}) && size(varargin{1},1) ~= size(varargin{2},1)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Dimension mismatch between center and generator matrix.'));                
                else
                    obj.Z = [varargin{1},varargin{2}]; 
                    obj.halfspace = [];
                end
            
            elseif nargin > 2
                
                % too many input arguments
                throw(CORAerror('CORA:tooManyInputArgs',2));
            end
        end
        
        % set parent object properties
        obj.dimension = size(obj.Z,1);
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random zonotope
    Z = enclosePoints(points,varargin) % enclose point cloud with zonotope
end

end

%------------- END OF CODE --------------