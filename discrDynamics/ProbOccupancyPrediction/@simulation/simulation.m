function Obj = simulation(varargin)
% simulation - ???
% Purpose:  1. Object constructor
%           2. Copy constructor
% Pre:      1st Parameter - simulation options
%           2nd Parameter - markov chain specification
%           Object as Parameter - Copy constructor
% Post:     Return a created object
%
% Syntax:
%    Obj = simulation(varargin)
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
% Last update:   01-October-2006
%                17-August-2007
%                17-June-2008
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% ------------------------------ END OF CODE ------------------------------
