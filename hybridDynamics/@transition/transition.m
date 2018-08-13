function Obj = transition(varargin)
% transition - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - guard: contSet
%    input2 - reset function (only linear map!): Ax+b, with struct:
%    reset.A, reset.b
%    input3 - target: int (id of target location)
%    input4 - input label: char array
%    input5 - output label: char array
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      02-May-2007 
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

% If no argument is passed (default constructor)
if nargin == 0
    disp('Transition needs more input values');
    Obj.guard=[];  % also add halfspace representation
    Obj.reset=[];
    Obj.target=[];
    Obj.inputLabel=[];
    Obj.outputLabel=[];
    % Register the variable as an object
    Obj = class(Obj, 'transition');    
    
% If 5 arguments are passed
elseif nargin == 5
    %List elements of the class
    %Obj.id=nextID('transition');  
    Obj.guard = varargin{1};
%     if isa(varargin{1},'halfspace') || isa(varargin{1},'constrainedHyperplane')
%         Obj.guard=varargin{1};
%     else
%         Obj.guard=polytope(varargin{1});  % also add halfspace representation
%     end
    Obj.reset=varargin{2};
    Obj.target=varargin{3};
    Obj.inputLabel=varargin{4};
    Obj.outputLabel=varargin{5};

    % Register the variable as an object
    Obj = class(Obj, 'transition');
        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'transition')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------