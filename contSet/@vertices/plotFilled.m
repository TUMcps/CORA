function plotFilled(varargin)
% plotFilled - Plots convex hull of vertices that are projected onto the first
% and second coordinate; plots are filled with color
%
% Syntax:  
%    plotFilled(V)
%
% Inputs:
%    V - vertices object 
%
% Outputs:
%    none
%
% Example: 
%    V=vertices(rand(2,6));
%    plotFilled(V)
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
    V=varargin{1};
    type='b';
    
%If two arguments are passed    
elseif nargin>=2
    V=varargin{1};
    type(1:length(varargin)-1)=varargin(2:end);
end

% convert vertices to polygon
p = polygon(V);

% plot
if ~isempty(p)
    fill(p(1,:),p(2,:),type{:});  
end


%------------- END OF CODE --------------