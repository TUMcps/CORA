function Obj = road(varargin)
% Purpose:  1. Object constructor
%           2. Copy constructor
% Pre:      1st Parameter - road width
%           2nd Parameter - segment length
%           3rd Parameter - discretization (angle, x-offset, y-offset)
% Post:     Return a created object
% Built:    16.11.06,MA
% Modified: 21.11.07,MA


% If no argument is passed (default constructor)
if nargin == 0
    disp('Road needs more input values');
    Obj=[];
    % Register the variable as an object
    Obj = class(Obj, 'road');        
    
% If two arguments are passed 
elseif nargin == 3
    %=======================================================
    Obj.width=varargin{1};
    Obj.segmentLength=varargin{2};
    Obj.nrOfDevSegments=varargin{3};          
    Obj.segments=[];
    Obj.prototypeVertices=[];
    %=======================================================
    % Register the variable as an object
    Obj = class(Obj, 'road');   
    
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'road')
    Obj = varargin{1};
    
% Otherwise use a specific constructor    
else
    disp('Road needs more/less input values');
    Obj=[];
end