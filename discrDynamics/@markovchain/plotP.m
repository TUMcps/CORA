function plotP(obj,varargin)
% plotP - plots the 2D probability distribution of a Markov chain
%
% Syntax:  
%    plotP(varargin)
%
% Inputs:
%    obj - markovchain object
%    varargin{1} - probability vector
%    varargin{2} - color
%    varargin{3} - plotStyle
%    varargin{4} - newfigure
%
% Outputs:
%    -
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      17-August-2007
% Last update:  14-November-2007 (implemented different colors)
%               26-March-2008 (implemented polytope plot)
% Last revision:---


%------------- BEGIN CODE --------------

partition = obj.field;

%no color specified
if nargin==2
    p=varargin{1};
    color='b'; %b=blue
    plotStyle='fill';
    newfigure=0;
    
%color specified
elseif nargin==3
    p=varargin{1};
    color=varargin{2};
    plotStyle='fill';
    newfigure=0;
    
%color and newfigure specified
elseif nargin==5
    p=varargin{1};
    color=varargin{2};
    plotStyle=varargin{3};
    newfigure=varargin{4};  
end

if newfigure
    figure; 
end

l=linspace(1,0,100)';
o=ones(100,1);

%prepare colormap
switch color
  case 'b' %blue
    colormap([l,l,o]);
  case 'r' %red
    colormap([o,l,l]);
  case 'g' %green
    colormap([l,o,l]);
  case 'k' %black
    colormap([l,l,l]);    
  otherwise
    disp('Error: unspecified color');
end


hold on
%get maximum probability for normalization
pMax=max(p);
%find nonzero probabilities
ind=find(p);
for i=1:length(ind)
    cellNr=ind(i)-1;
    if cellNr~=0
        %generate polytope out of cell
        IHP=cellPolytopes(partition,cellNr); %IHP:interval hull polytope        
        if strcmp(plotStyle,'polytope')
            options.color=color;
            options.linestyle='none';
            options.shade=p(cellNr+1)/pMax;
            plot(IHP{1},options);            
        else
            %get vertices
            V=vertices(IHP{1});
            x=V(1,:);
            y=V(2,:);
            k=convhull(x,y);
            fill(x(k),y(k),0,'EdgeColor','none'); %workaround to adjust colorbar
            fill(x(k),y(k),p(cellNr+1),'EdgeColor','none');
        end
    end
end

%------------- END OF CODE --------------