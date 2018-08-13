function Obj = probZonotope(varargin)
% probZonotope - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = zonotope(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - zonotope matrix
%    input2 - probabilistic generators
%    input3 - weighting, mean and variance row vectors
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff
% Written:      03-August-2007 
% Last update:  26-February-2008
%               20-March-2015
% Last revision: ---

%------------- BEGIN CODE --------------

%Probabilistic Zonotope is superior to class zonotope
superiorto('zonotope');
superiorto('interval');

% If no argument is passed
if nargin == 0
    disp('This class needs more input values');
    Obj=[];
       
    
% If 3 arguments are passed
elseif nargin == 3
    %List elements of the class
    Obj.Z=varargin{1}; 
    Obj.g=varargin{2}; 
    Obj.cov=[]; %covariance matrix
    Obj.gauss=0; %flag, determining if Obj.cov is updated
    Obj.gamma=varargin{3}; %cut-off mSigma value
    %Generate parent object
    cSet=contSet(length(varargin{1}(:,1)));        

    % Register the variable as an child object of contSet
    Obj = class(Obj, 'probZonotope', cSet); 
    
    %Update covariance matrix
    Obj.cov=sigma(Obj);
    Obj.gauss=1; %flag, determining if Obj.cov is updated
    
    
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'probZonotope')
    Obj = varargin{1};
    
% Else if not enough or too many inputs are passed    
else
    disp('This class needs more/less input values');
    Obj=[];
end

%------------- END OF CODE --------------