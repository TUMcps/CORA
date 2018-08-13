function obj = mptPolytope(varargin)
% mptPolytope - object and copy constructor 
%
% Syntax:  
%    object constructor: obj = zonotope(varargin)
%    copy constructor: obj = otherObj
%
% Inputs:
%    V - vertices
%    C/d - halfspace representation
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalhull,  polytope

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

%Zonotope is superior to class interval
superiorto('interval','intervalMatrix','matZonotope');


% If no argument is passed (default constructor)
if nargin == 0
    disp('Zonotope needs more input values');
    obj.P=[];
    obj.halfspace=[];
    
    % Register the variable
    obj = class(obj, 'mptPolytope');

    
% If 1 argument is passed
elseif nargin == 1
    %List elements of the class
    try %MPT V3
        obj.P=Polyhedron(varargin{1});
    catch %MPT V2
        obj.P=polytope(varargin{1});
    end
    obj.halfspace=[];

    % Register the variable as an child object of contSet
    obj = class(obj, 'mptPolytope'); 
    
% If 2 arguments are passed
elseif nargin == 2
    %List elements of the class
    try %MPT V3
        obj.P=Polyhedron(varargin{1},varargin{2});
    catch %MPT V2
        obj.P=polytope(varargin{1},varargin{2});
    end
    obj.halfspace=[];

    % Register the variable as an child object of contSet
    obj = class(obj, 'mptPolytope'); 
    
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'mptPolytope')
    obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    obj=[];
end

%------------- END OF CODE --------------