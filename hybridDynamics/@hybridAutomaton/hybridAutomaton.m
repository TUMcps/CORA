function Obj = hybridAutomaton(varargin)
% hybridAutomaton - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - location array
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
% Written: 03-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


% If no argument is passed (default constructor)
if nargin == 0
    disp('HybridAutomaton needs more input values');
    Obj=[];
    % Register the variable as an object
    Obj = class(Obj, 'hybridAutomaton');    
    
% If 1 argument is passed
elseif nargin == 1
    %List elements of the class
    %Obj.id=nextID('hybridAutomaton');  
    Obj.location=varargin{1};
    Obj.result=[];

    % Register the variable as an object
    Obj = class(Obj, 'hybridAutomaton');
        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'hybridAutomaton')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------