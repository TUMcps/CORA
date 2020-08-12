classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonotope
% zonotope - Object and Copy Constructor 
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
%    zono = zonotope(c,G);
%    plot(zono,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      14-September-2006 
% Last update:  22-March-2007
%               04-June-2010
%               08-February-2011
%               18-November-2015
%               05-December-2017 (DG) class is redefined in complience with
%               the new standard.
%               28-April-2019 code shortened
%               1-May-2020 (NK) new constructor + removed orientation prop.
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    Z (:,:) {mustBeNumeric} = []; % zonotope center and generator Z = [c,g_1,...,g_p]
    halfspace = [];     % halfspace representation of the zonotope
    contSet = [];
end

methods

    function Obj = zonotope(varargin)
        
        % If no argument is passed (default constructor)
        if nargin == 0
            Obj.Z = [];
            Obj.halfspace = [];

            % Generate parent object
            Obj.contSet = contSet();

        % If 1 argument is passed
        else
            
            % Generate parent object
            if ~isempty(varargin{1})
                Obj.contSet = contSet(length(varargin{1}(:,1)));
            else
                Obj.contSet = contSet();
            end
            
            if nargin == 1
            
                % input is a zonotope -> copy object
                if isa(varargin{1},'zonotope')
                    Obj = varargin{1};
                else
                    % List elements of the class
                    Obj.Z = varargin{1}; 
                    Obj.halfspace = [];

                end

            % If 2 arguments are passed
            elseif nargin == 2

                % List elements of the class
                Obj.Z = [varargin{1},varargin{2}]; 
                Obj.halfspace = [];
            
            end
        
        end
        
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random zonotope
    Z = enclosePoints(points,varargin) % enclose point cloud with zonotope
end


end
%------------- END OF CODE --------------