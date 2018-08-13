function Obj = contSet(varargin)
% contSet - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - dimension: int
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
    disp('ContSet needs more input values');
    Obj.id=0;
    Obj.dimension=0;
    % Register the variable as an object
    Obj = class(Obj, 'contSet');    
    
% If 1 argument is passed
elseif nargin == 1
    %List elements of the class
    %Obj.id=nextID('contSet');  
    Obj.id=0;
    Obj.dimension=varargin{1};

    % Register the variable as an object
    Obj = class(Obj, 'contSet');
        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'contSet')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------