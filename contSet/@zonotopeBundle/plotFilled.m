function plotFilled(varargin)
% plotFilled - Plots 2-dimensional over-approximative projection of a 
% zonotope bundle and fills it out
%
% Syntax:  
%    plot(varargin)
%
% Inputs:
%    Z - zonotope bundle object
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

% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    Zbundle=varargin{1};
    dimensions=[1,2];
    type='b';
    
%If two arguments are passed    
elseif nargin==2
    Zbundle=varargin{1};
    dimensions=varargin{2};
    type='b';
   
%If three or more arguments are passed
elseif nargin>=3
    Zbundle=varargin{1};
    dimensions=varargin{2};   
    type(1:length(varargin)-2)=varargin(3:end);
end


%compute polytopes
for i=1:Zbundle.parallelSets
    %delete zero generators
    Z=deleteZeros(Zbundle.Z{i});
    %project zonotope
    Z=project(Z,dimensions);
    %convert to polytope
    P{i}=polytope(Z);
end

%intersect polytopes
Pint=P{1};
for i=2:Zbundle.parallelSets
    Pint=Pint&P{i};
end
    
%Plot polytope
plotFilled(Pint,[1 2],type{:});

%------------- END OF CODE --------------