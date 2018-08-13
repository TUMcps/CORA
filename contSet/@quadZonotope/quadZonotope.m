function Obj = quadZonotope(varargin)
% quadZonotope - Object and Copy Constructor 
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
% Written:      04-September-2012
% Last update:  01-January-2016
%               26-January-2016
% Last revision:---

%------------- BEGIN CODE --------------

%Zonotope is superior to class interval
superiorto('interval','intervalMatrix','matZonotope');


% If no argument is passed (default constructor)
if nargin == 0
    disp('quadZonotope needs more input values');
    Obj.Z=[];
    Obj.O=[]; %orientation
    Obj.halfspace=[];
    
    %Generate parent object
    cSet=contSet(0);
    
    % Register the variable as an child object of contSet
    Obj = class(Obj, 'quadZonotope', cSet);

% If 1 argument is passed
elseif nargin == 1
    %List elements of the class
    Obj.c=varargin{1};
    Obj.G=[]; 
    Obj.Gsquare=[]; 
    Obj.Gquad=[];
    Obj.Grest=[]; 
    Obj.O=[]; 
    Obj.halfspace=[];
    
    %Generate parent object
    if ~isempty(varargin{1})
        cSet=contSet(length(varargin{1}(:,1)));
    else
        cSet=contSet(0);
    end

    % Register the variable as an child object of contSet
    Obj = class(Obj, 'quadZonotope', cSet);     
    
% If 5 arguments are passed
elseif nargin == 5
    %List elements of the class
    Obj.c=varargin{1};
    Obj.G=varargin{2}; 
    if isempty(varargin{3})
        Obj.Gsquare=zeros(length(Obj.c), length(Obj.G(1,:)));
    else
        Obj.Gsquare=varargin{3}; 
    end
    if isempty(varargin{4})
        if length(Obj.G(1,:))>1
            % binom(length(Obj.G(1,:)),2)
            n = length(Obj.G(1,:));
            gens = n*(n-1)/2; 
            Obj.Gquad=zeros(length(Obj.c), int64(gens));
        else
            Obj.Gquad=[];
        end
    else
        Obj.Gquad=varargin{4}; 
    end
    if isempty(varargin{5})
        Obj.Grest=zeros(length(Obj.c), 1); 
    else
        Obj.Grest=varargin{5}; 
    end
    Obj.O=[]; 
    Obj.halfspace=[];
    
    %Generate parent object
    if ~isempty(varargin{1})
        cSet=contSet(length(varargin{1}(:,1)));
    else
        cSet=contSet(0);
    end

    % Register the variable as an child object of contSet
    Obj = class(Obj, 'quadZonotope', cSet); 
    
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'quadZonotope')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------
