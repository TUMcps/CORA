function Obj = binds( varargin )
%COMPONENT - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - input id
%    input2 - component id
%    input3 - state variable id
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Victor CHARLENT
% Written:      24-May-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% If no argument is passed
if nargin == 0
    disp('This class needs more input values');
    Obj.input=[];
    Obj.component=[];
    Obj.state=[];
    
    % Register the variable as an object
    Obj = class(Obj, 'bind');
    
% If 3 arguments are passed
elseif nargin == 3
    %List elements of the class
    Obj.input=varargin{1};
    Obj.component = varargin{2};
    Obj.state=varargin{3};

    % Register the variable as an object
    Obj = class(Obj, 'bind');
        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'bind')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

end



