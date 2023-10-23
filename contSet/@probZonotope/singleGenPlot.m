function singleGenPlot(probZ,varargin)
% singleGenPlot - Plots 2-dimensional projection of a probabilistic
%    zonotope with a maximum of 5 generators
%
% Syntax:
%    singleGenPlot(probZ,dimensions)
%
% Inputs:
%    probZ - probabilistic zonotope object
%    dimensions - dimensions that should be projected (optional) 
%
% Outputs:
%    ---
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    singleGenPlot(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       03-August-2007
% Last update:   29-February-2008
%                02-September-2009
%                04-September-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin <= 3
    % parse input arguments
    [type,m] = setDefaultValues({'solid',probZ.gamma},varargin);

    % check input arguments
    inputArgsCheck({{probZ,'att','probZonotope'};
                    {type,'str',{'solid','mesh'}};
                    {m,'att','numeric','nonnan'}});  
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

%dimension of the single generator
n = dim(probZ);

%init number of plotted points
nrOfPoints=1e3;

%get center
c = center(probZ);

if n==1
    %compute Sigma
    Sigma=sigma(probZ);    
    
    if length(probZ.Z)==1
        x=linspace(-m*norm(probZ.g),m*norm(probZ.g),nrOfPoints);
        for i=1:nrOfPoints    
            f(i)=gaussian(x(i),Sigma);
        end
    else
        c1=-sum(abs(probZ.Z(2:end)));
        c2=-c1;
        l1=linspace(-m*norm(probZ.g),c1,nrOfPoints);
        l2=linspace(c2,m*norm(probZ.g),nrOfPoints);
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
    Sigma=norm(probZ.g)^2;  
    l=linspace(-m,m,nrOfPoints);
    for i=1:nrOfPoints
        x(i)=c(1)+probZ.g(1,1)*l(i);
        y(i)=c(2)+probZ.g(2,1)*l(i);
        f(i)=gaussian(norm(probZ.g)*l(i),Sigma);
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
    norm(probZ.g)*(l(2)-l(1))*sum(f)
end

% ------------------------------ END OF CODE ------------------------------
