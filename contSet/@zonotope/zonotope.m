classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonotope
% zonotope - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = zonotope(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - zonotope matrix
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalhull,  polytope

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
%               04-June-2010
%               08-February-2011
%               18-November-2015
%               05-December-2017 (DG) class is redefined in complience with
%               the new standard.
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    Z=[];
    O=[];
    halfspace=[];
    contSet = [];
end
   
methods

    function Obj = zonotope(varargin)
        
        %Obj.contSet = Obj@contSet(length(varargin{1}(:,1)));
        % If no argument is passed (default constructor)
        if nargin == 0
            disp('Zonotope needs more input values');
            Obj.Z=[];
            Obj.O=[]; %orientation
            Obj.halfspace=[];

            %Generate parent object
            %cSet=contSet(0); DEL !!!
            Obj.contSet = contSet(0);

            % Register the variable as an child object of contSet
            %Obj = class(Obj, 'zonotope', cSet); % !!!!!


        % If 1 argument is passed
        elseif nargin == 1
            %is input a zonotope?
            if isa(varargin{1},'zonotope')
                Obj = varargin{1};
            else
                %List elements of the class
                Obj.Z=varargin{1}; 
                Obj.O=[]; 
                Obj.halfspace=[];

                %Generate parent object
                if ~isempty(varargin{1})
                    %cSet=contSet(length(varargin{1}(:,1)));
                    Obj.contSet = contSet(length(varargin{1}(:,1)));
                else
                    %cSet=contSet(0);
                    Obj.contSet = contSet(0);
                end

                % Register the variable as an child object of contSet
                
                %Obj = class(Obj, 'zonotope', cSet); !!!!!!
            end


        % If 2 arguments are passed
        elseif nargin == 2
            %List elements of the class
            Obj.Z=varargin{1}; 
            Obj.O=varargin{2};  
            Obj.halfspace=[];

            %Generate parent object
            if ~isempty(varargin{1})
                Obj.contSet=contSet(length(varargin{1}(:,1)));
            else
                Obj.contSet=contSet(0);
            end

            % Register the variable as an child object of contSet
            %Obj = class(Obj, 'zonotope', cSet);  

        % Else if the parameter is an identical object, copy object    
        elseif isa(varargin{1}, 'zonotope')
            Obj = varargin{1};

        % Else if not enough or too many inputs are passed    
        else
            disp('This class needs more/less input values');
            Obj=[];
        end
    end
end
end
%------------- END OF CODE --------------