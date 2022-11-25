function Obj = markovchain(varargin)
% Purpose:  1. Object constructor
%           2. Copy constructor
% Pre:      1st Parameter - partition
%           2nd Parameter - transition matrix
%           Object as Parameter - Copy constructor
% Post:     Return a created object
% Tested:   14.09.06,MA
% Modified: 17.08.07,MA

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

%-----------------------------------------------------------
function [T]=init_T(field)
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