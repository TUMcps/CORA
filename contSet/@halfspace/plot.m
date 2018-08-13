function plot(varargin)
% plot - Plots 2-dimensional projection of a halfspace
%
% Syntax:  
%    plot(obj,dimensions,type)
%
% Inputs:
%    obj - halfspace object
%    dimensions - dimensions that should be projected (optional); assume
%    that other entries of the normal vector are zeros
%    type - plot type
%    
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

% Author:       Matthias Althoff
% Written:      23-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    obj=varargin{1};
    dimensions=[1,2];
    type='frame';
    
%If two arguments are passed    
elseif nargin==2
    obj=varargin{1};
    dimensions=varargin{2};
    type='frame';
    
%If too many arguments are passed
elseif nargin==3
    obj=varargin{1};
    dimensions=varargin{2};   
    type=varargin{3};
end

% compute slope
m = -obj.c(dimensions(1))/obj.c(dimensions(2));
b = obj.d/obj.c(dimensions(2));

%plot reference line
refline(m,b);

%------------- END OF CODE --------------