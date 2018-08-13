function singleGenPlot(varargin)
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

% Author:       Matthias Althoff
% Written:      03-August-2007
% Last update:  29-February-2008
%               02-September-2009
%               04-September-2009
% Last revision: ---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    pZ=varargin{1};
    type='solid';
    m=pZ.gamma;
    
%If two arguments are passed    
elseif nargin==2
    pZ=varargin{1};
    type=varargin{2}; 
    m=pZ.gamma;
    
%If three arguments are passed    
elseif nargin==3
    pZ=varargin{1};
    type=varargin{2}; 
    m=varargin{3};     
    
%If too many arguments are passed
else
    disp('Error: too many inputs');
    pZ=varargin{1};
    type=varargin{2};    
end

%dimension of the single generator
dim=length(pZ.g);

%init number of plotted points
nrOfPoints=1e3;

%get center
c=pZ.Z(:,1);

if dim==1
    %compute Sigma
    Sigma=sigma(pZ);    
    
    if length(pZ.Z)==1
        x=linspace(-m*norm(pZ.g),m*norm(pZ.g),nrOfPoints);
        for i=1:nrOfPoints    
            f(i)=gaussian(x(i),Sigma);
        end
    else
        c1=-sum(abs(pZ.Z(2:end)));
        c2=-c1;
        l1=linspace(-m*norm(pZ.g),c1,nrOfPoints);
        l2=linspace(c2,m*norm(pZ.g),nrOfPoints);
        for i=1:nrOfPoints    
            x(i)=l1(i);
            f(i)=gaussian(l1(i)-c1,Sigma);
        end    
        for i=(nrOfPoints+1):2*nrOfPoints   
            x(i)=l2(i-nrOfPoints);
            f(i)=gaussian(l2(i-nrOfPoints)-c2,Sigma);
        end         
    end
    plot(c+x,f);
    (x(2)-x(1))*sum(f)
else
    Sigma=norm(pZ.g)^2;  
    l=linspace(-m,m,nrOfPoints);
    for i=1:nrOfPoints
        x(i)=c(1)+pZ.g(1,1)*l(i);
        y(i)=c(2)+pZ.g(2,1)*l(i);
        f(i)=gaussian(norm(pZ.g)*l(i),Sigma);
    end

    xMin=min(x);
    xMax=max(x);
    yMin=min(y);
    yMax=max(y);

    xRange=linspace(xMin,xMax,20);
    yRange=linspace(yMin,yMax,20);

    colormap([0,0,1]);

    %plot graph
    if strcmp(type,'mesh')
        mesh(xRange,yRange,zeros(20));
        hold on
    %     surf(x,y,diag(f),'FaceColor','interp',...
    %     'EdgeColor','none',...
    %     'FaceLighting','phong')    
        plot3(x,y,f)
        %hidden off
    else
        surf(x,y,diag(f),'FaceColor','interp',...
        'EdgeColor','none',...
        'FaceLighting','phong')
        material dull
        %alpha(.4)
        %set lights
        %camlight right
        %camlight left
        camlight headlight
    end
    norm(pZ.g)*(l(2)-l(1))*sum(f)
end


    
    

%------------- END OF CODE --------------