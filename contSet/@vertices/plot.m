function plot(varargin)
% plot - Plots convex hull of vertices that are projected onto the first
% and second coordinate
%
% Syntax:  
%    plot(V)
%
% Inputs:
%    V - vertices object 
%
% Outputs:
%    none
%
% Example: 
%    V=vertices(rand(2,6));
%    plot(V)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      30-September-2006 
% Last update:  22-March-2007
%               18-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    V=varargin{1};
    type='frame';
    
%If two arguments are passed    
elseif nargin==2
    V=varargin{1};
    type=varargin{2};
    
%If three arguments are passed    
elseif nargin==3
    V=varargin{1};
    type=varargin{2};  
    val=varargin{3};  
    
%If too many arguments are passed
else
    disp('Error: too many inputs');
    V=varargin{1};
    type=varargin{2};    
end

% convert vertices to polygon
p = polygon(V);

if ~isempty(p)
    hold on
    %try to use linespec
    try
        plot(p(1,:),p(2,:),type);  
    %use one of the special definitions
    catch
        if strcmp(type,'frame')
            plot(p(1,:),p(2,:),'b-');  
        elseif strcmp(type,'filledFrame')
            fill(p(1,:),p(2,:),'w','EdgeColor','b'); 
        elseif strcmp(type,'greenFrame')
            fill(p(1,:),p(2,:),'w','EdgeColor','g');            
        elseif strcmp(type,'lightgray')
            fill(p(1,:),p(2,:),[.75 .75 .75],'EdgeColor','none');
        elseif strcmp(type,'gray')
            fill(p(1,:),p(2,:),[.675 .675 .675],'EdgeColor','none');    
        elseif strcmp(type,'darkgray')
            fill(p(1,:),p(2,:),[.6 .6 .6],'EdgeColor','none');
        elseif strcmp(type,'black')
            fill(p(1,:),p(2,:),[.0 .0 .0],'EdgeColor','none');               
        elseif strcmp(type,'grayEdge')
            plot(p(1,:),p(2,:),'Color',[.8 .8 .8]); 
        elseif strcmp(type,'grayEdge2')
            plot(p(1,:),p(2,:),'Color',[.8 .8 .8],'lineWidth',1); 
        elseif strcmp(type,'grayEdgeThick')
            plot(p(1,:),p(2,:),'Color',[.8 .8 .8],'lineWidth',3);
        elseif strcmp(type,'blackEdge')
            plot(p(1,:),p(2,:),'k-','lineWidth',2);   
        elseif strcmp(type,'blackEdgeThick')
            plot(p(1,:),p(2,:),'k-','lineWidth',3);   
        elseif strcmp(type,'whiteEdge')
            plot(p(1,:),p(2,:),'w-','lineWidth',2);
        elseif strcmp(type,'whiteEdgeDashed')
            plot(p(1,:),p(2,:),'w--','lineWidth',2);
        elseif strcmp(type,'blackFrame')
            fill(p(1,:),p(2,:),'w','EdgeColor','k'); 
            %plot(p(1,:),p(2,:),'k-');           
        elseif strcmp(type,'grayFrame')
            fill(p(1,:),p(2,:),[.8 .8 .8],'EdgeColor','k');  
        elseif strcmp(type,'grayTones')
            l=linspace(1,0,100)';
            colormap([l,l,l]); 
            fill(p(1,:),p(2,:),val,'EdgeColor','none');            
        elseif strcmp(type,'vertices')
            plot(p(1,:),p(2,:),'k+');   
        elseif strcmp(type,'2dIn3d')
            zCoordinates=10*ones(length(p(1,:)),1);
            plot3(p(1,:),p(2,:),zCoordinates,'k-');             
        end
    end
end


%------------- END OF CODE --------------