function [GM]=plot(varargin)
% plot - Plots 2-dimensional projection of a zonotope with a maximum of 5
% generators
%
% Syntax:  
%    plot(Z,dimensions)
%
% Inputs:
%    Z - zonotope object
%    dimensions - dimensions that should be projected (optional) 
%
% Outputs:
%    none
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 03-August-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    pZ=varargin{1};
    type='solid';
    dimensions=[1,2];
    m=pZ.gamma;
    
%If two arguments are passed    
elseif nargin==2
    pZ=varargin{1};
    type=varargin{2};
    dimensions=[1,2];
    m=pZ.gamma;
    
%If three arguments are passed    
elseif nargin==3
    pZ=varargin{1};
    type=varargin{2};
    dimensions=varargin{3};
    m=pZ.gamma;
    
%If four arguments are passed    
elseif nargin==4
    pZ=varargin{1};
    type=varargin{2};
    dimensions=varargin{3};
    m=varargin{4};     
    
%If too many arguments are passed
else
    disp('Error: too many inputs');
    pZ=varargin{1};
    dimensions=varargin{2};    
end
    

%compute enclosing probability
[eP] = enclosingProbability(pZ,m,dimensions);

%open new figure
%figure
%colormap('cool')
%colormap('hot')

% l=linspace(1,0,100)';
% o=ones(100,1);
% colormap([l,l,o]);

%colormap([0,0,0]);


%plot graph
if strcmp(type,'mesh')
    mesh(eP.X,eP.Y,eP.P);
    hidden off
elseif strcmp(type,'meshHide')
    colormap([0,0,0]);
    mesh(eP.X,eP.Y,eP.P); 
elseif strcmp(type,'blue')
    colormap([0,0,1]);
    surf(eP.X,eP.Y,eP.P,'FaceColor','b',...
    'EdgeColor','none',...
    'FaceLighting','phong')
    material dull
    %alpha(.4)
    %set lights
    %camlight right
    %camlight left
    camlight headlight   
elseif strcmp(type,'green')
    surf(eP.X,eP.Y,eP.P,'FaceColor','g',...
    'EdgeColor','none',...
    'FaceLighting','phong')
    material dull
    %alpha(.4)
    %set lights
    %camlight right
    %camlight left
    camlight headlight    
elseif strcmp(type,'dark')
    surf(eP.X,eP.Y,eP.P,'FaceColor',[0.2 0.2 0.2],...
    'EdgeColor','none',...
    'FaceLighting','phong')
    material dull
    %alpha(.4)
    %set lights
    %camlight right
    %camlight left
    camlight headlight   
elseif strcmp(type,'light')
    surf(eP.X,eP.Y,eP.P,'FaceColor',[0.5 0.5 0.5],...
    'EdgeColor','none',...
    'FaceLighting','phong')
    material dull
    %alpha(.4)
    %set lights
    %camlight right
    %camlight left
    camlight headlight        
else
    l=linspace(1,0,100)';
    o=ones(100,1);
    colormap([l,l,l]);

    cmap = colormap;
    newCmap=[ones(1,3);cmap(2:end,:)];
    %newCmap=cmap;
    colormap(newCmap);
    
    surf(eP.X,eP.Y,eP.P,'FaceColor','interp',...
    'EdgeColor','none')    
    
%     surf(eP.X,eP.Y,eP.P,'FaceColor','interp',...
%     'EdgeColor','none',...
%     'FaceLighting','phong')
%     material dull 
end



%daspect([40 40 1])
%axis([-2 1.8 -2 1 0 1.3])
% xlim([-2 1.8])
% ylim([-2 1])
% zlim([0 1.3])

%------------- END OF CODE --------------