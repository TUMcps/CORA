function Obj = component( varargin )
%COMPONENT - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - name of the component
%    input2 - collection of components
%    input3 - collection of locations
%    input4 - collection of bind
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
    Obj.name=[];
    Obj.components = [];
    Obj.loc=[];
    Obj.binds=[];
    
    % Register the variable as an object
    Obj = class(Obj, 'component');
    
% If 5 arguments are passed
elseif nargin == 4
    %List elements of the class
    Obj.name=varargin{1};
    Obj.components = varargin{2};
    Obj.loc=varargin{3};
    Obj.binds = varargin{4};

    % Register the variable as an object
    Obj = class(Obj, 'component');

        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'component')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end






end

