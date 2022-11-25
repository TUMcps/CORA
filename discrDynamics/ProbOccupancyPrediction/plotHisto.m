function plotHisto(varargin)
% plotHisto - plots the histogram of the position probability function
%
% Syntax:  
%    plotHisto(obj,posProb)
%
% Inputs:
%    obj - partition object
%    posProb - position probability vector
%
% Outputs:
%    ---
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      22-August-2008
% Last update:  19-October-2009
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==2
    obj=varargin{1};
    posProb=varargin{2};
    plotStyle='-k';
elseif nargin==3
    obj=varargin{1};
    posProb=varargin{2};
    plotStyle=varargin{3};    
else
    disp('more inputs needed');
end

%get position interval length
posInt=obj.nrOfSegments(1);

%generate position vector
posVec=linspace(obj.intervals(1,1),obj.intervals(1,2),posInt+1);

%In order to plot the probability DENSITY function, one has to devide the
%probabilities by the position interval length
%posProb=posProb/(posVec(2)-posVec(1));


%figure;

%init x,y
x(1)=posVec(1);
y(1)=posProb(1);
x(2)=posVec(2);
y(2)=y(1);

%loop
for i=1:obj.nrOfSegments(1)
    x(2*i-1)=posVec(i);
    y(2*i-1)=posProb(i);
    x(2*i)=posVec(i+1);
    y(2*i)=y(2*i-1);
end

%plot
plot(x,y,plotStyle);

%------------- END OF CODE --------------