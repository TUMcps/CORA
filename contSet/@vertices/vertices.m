function Obj = vertices(varargin)
% vertices - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = vertices(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - Matrix with points as column vectors
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope,  intervalhull

% Author: Matthias Althoff
% Written: 22-March-2007 
% Last update: 22-March-2007
% Last revision: 22-March-2007

%------------- BEGIN CODE --------------

% If no argument is passed (default constructor)
if nargin == 0
    disp('Vertices needs more input values');
    Obj=[];
    % Register the variable as an object
    Obj = class(Obj, 'vertices');    
    
% If n arguments are passed
elseif nargin == 1
    %List elements of the class
    Obj.V=varargin{1}; 

    % Register the variable as an object
    Obj = class(Obj, 'vertices');  
    
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'vertices')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------
