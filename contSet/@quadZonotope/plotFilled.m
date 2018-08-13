function plotFilled(varargin)
% plotFilled - Plots 2-dimensional overapproximation of a quadZonotope
%
% Syntax:  
%    plot(Z,dimensions,type,splits)
%
% Inputs:
%    qZ - quadZonotope object
%    dimensions - dimensions that should be projected (optional) 
%    type - plot type
%    splits - number of splits for refinement
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
% Written:      12-August-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    qZ = varargin{1};
    dimensions = [1,2];
    splits = 0;
    order = [];
    type = 'b';
    
%If two arguments are passed    
elseif nargin==2
    qZ = varargin{1};
    dimensions = varargin{2};
    splits = 0;
    order = [];
    type = 'b';
    
%If three arguments are passed
elseif nargin==3
    qZ = varargin{1};
    dimensions = varargin{2};   
    splits = varargin{3};
    order = [];
    type = 'b';
    
%If four arguments are passed
elseif nargin==4
    qZ = varargin{1};
    dimensions = varargin{2};   
    splits = varargin{3}; 
    order = varargin{4};
    type = 'b';
      
%If five arguments are passed
elseif nargin>=5
    qZ = varargin{1};
    dimensions = varargin{2};   
    splits = varargin{3};  
    order = varargin{4};
    type(1:length(varargin)-4) = varargin(5:end);
end


%initilize splitting 
qZsplit{1} = qZ;

for i=1:splits
    qZnew = [];
    for j=1:length(qZsplit)
        res = splitLongestGen(qZsplit{j});
        qZnew{end+1} = res{1};
        qZnew{end+1} = res{2};
    end
    qZsplit = qZnew;
end


%plot
for i=1:length(qZsplit)
    qZproj = project(qZsplit{i}, dimensions);
    Zproj = zonotope(qZproj);
    
    %reduction
    if ~isempty(order)
        Zproj = reduce(Zproj,'girard',order);
    end
    plotFilled(Zproj,[1 2],type{:});
end

%------------- END OF CODE --------------