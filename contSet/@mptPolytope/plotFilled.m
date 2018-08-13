function plotFilled(varargin)
% plotFilled - Plots 2-dimensional projection of a mptPolytope and fills it
% with a color
%
% Syntax:  
%    plotFilled(P,dimensions,type)
%
% Inputs:
%    P - mptPolytope
%    dimensions - dimensions that should be projected (optional) 
%    type - plot type (optional) 
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
% Written:      18-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    P=varargin{1};
    dimensions=[1,2];
    type='b';
    
%If two arguments are passed    
elseif nargin==2
    P=varargin{1};
    dimensions=varargin{2};
    type='b';
    
%If too many arguments are passed
elseif nargin>=3
    P=varargin{1};
    dimensions=varargin{2};   
    type(1:length(varargin)-2)=varargin(3:end);
end
 
%compute vertices
V=vertices(P);
%Plot convex hull of projected vertices
if iscell(V)
    for i=1:length(V)
        plotFilled(vertices(V{i}(dimensions,:)),type{:});
    end
else
    try 
        plotFilled(vertices(V(dimensions,:)),type{:});
    catch
        disp('plot error');
    end
end


%------------- END OF CODE --------------