function plot(varargin)
% plot - Plots 2-dimensional projection of the zonotopes of the reachable
% set
%
% Syntax:  
%    plot(obj,dimensions)
%
% Inputs:
%    obj - linear interval system object
%    dimensions - dimensions that should be projected (optional) 
%
% Outputs:
%    none
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 30-April-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    obj=varargin{1};
    dimensions=[1,2];
    
%If two arguments are passed    
elseif nargin==2
    obj=varargin{1};
    dimensions=varargin{2};
    
%If too many arguments are passed
else
    disp('Error: too many inputs');
    obj=varargin{1};
    dimensions=varargin{2};    
end

%set maximum order of zonotopes that are plotted (otherwise, plot can be 
%extremly time consuming)
switch obj.dim
  case 2
    maxOrder=5;
  case 3
    maxOrder=4;
  case 4
    maxOrder=3;
  case 5
    maxOrder=2;    
  otherwise
    maxOrder=1; 
end
    
%plot each zonotope
for i=1:length(obj.reachSet)
    reachSet=reduce(obj.reachSet{i},'girard',maxOrder);
    plot(reachSet);
    hold on
end

%------------- END OF CODE --------------