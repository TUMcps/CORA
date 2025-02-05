function Obj = markovchain(varargin)
% markovchain - ???
% Purpose:  1. Object constructor
%           2. Copy constructor
% Pre:      1st Parameter - partition
%           2nd Parameter - transition matrix
%           Object as Parameter - Copy constructor
% Post:     Return a created object
%
% Syntax:
%    Obj = markovchain(varargin)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       14-September-2006
% Last update:   17-August-2007
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% If no argument is passed (default constructor)
if nargin == 0
    disp('Markovchain needs more input values');
    Obj=[];
    % Register the variable as an object
    Obj = class(Obj, 'markovchain');         
% If one argument is passed
elseif nargin == 1
    %======================================================
    Obj.field=varargin{1}; 
    Obj.T=[];
    %======================================================
    % Register the variable as an object
    Obj = class(Obj, 'markovchain');    
% If two arguments are passed
elseif nargin == 2
    %======================================================
    Obj.field=varargin{1}; 
    Obj.T=varargin{2};
    %======================================================
    % Register the variable as an object
    Obj = class(Obj, 'markovchain');        
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'markovchain')
    Obj = varargin{1};
% Otherwise use a specific constructor    
else
    disp('Markovchain needs more/less input values');
    Obj=[];
end


% Auxiliary functions -----------------------------------------------------

function [T]=aux_init_T(field)
% Purpose:  initializetransition matrices
% Pre:      field
% Post:     T


%initialize
nrOfSegments=get(field,'nrOfSegments');
totalNrOfSegments=prod(nrOfSegments);

T_init=zeros(totalNrOfSegments+1); %zero matrix
T_init(1,1)=1; %no transitions from outside area

T.T=T_init;
T.OT=T_init;

% ------------------------------ END OF CODE ------------------------------
