function Obj = simulation(varargin)
% Purpose:  1. Object constructor
%           2. Copy constructor
% Pre:      1st Parameter - simulation options
%           2nd Parameter - markov chain specification
%           Object as Parameter - Copy constructor
% Post:     Return a created object
% Tested:   14.09.06,MA
% Modified: 01.10.06,MA
% Modified: 17.08.07,MA
% Modified: 17.06.08,MA


% If no argument is passed (default constructor)
if nargin == 0
    disp('Simulation needs more input values');
    Obj=[];
    % Register the variable as an object
    Obj = class(Obj, 'simulation');    
    
% If 3 arguments are passed
elseif nargin == 2
    %======================================================
    Obj.simOptions=varargin{1};
    Obj.markovChainSpec=varargin{2};
    Obj.result=[];
    %======================================================
    % Register the variable as an object
    Obj = class(Obj, 'simulation');    
% Else if the parameter is an identical object, copy object    
elseif isa(varargin{1}, 'simulation')
    Obj = varargin{1};
% Otherwise use a specific constructor    
else
    disp('Simulation needs more/less input values');
    Obj=[];
end