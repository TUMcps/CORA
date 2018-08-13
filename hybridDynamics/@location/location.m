function Obj = location(varargin)
% location - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - name of the location: char array
%    input2 - invariant set: contSet
%    input3 - transition array: transition
%    input4 - continuous dynamics: contDynamics 
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 02-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

% If no argument is passed (default constructor)
if nargin == 0
    disp('Location needs more input values');
    Obj.name=[];
    Obj.id=[];
    Obj.invariant=[];  % also add halfspace representation
    Obj.transition=[];
    Obj.contDynamics=[];
    Obj.other=[];
    % Register the variable as an object
    Obj = class(Obj, 'location');    
    
% If 5 arguments are passed
elseif nargin == 5
    %List elements of the class
    %Obj.id=nextID('location');  
    Obj.name=varargin{1};
    Obj.id=varargin{2};
    if ~isempty(varargin{3})
        Obj.invariant=polytope(varargin{3});  % also add halfspace representation
    else
        Obj.invariant=[];
    end
    Obj.transition=varargin{4};
    Obj.contDynamics=varargin{5};
    Obj.other=[];

    % Register the variable as an object
    Obj = class(Obj, 'location'); 
        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'location')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------