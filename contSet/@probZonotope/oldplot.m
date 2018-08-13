function [xy,probXY]=oldplot(varargin)
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
    dimensions=[1,2];
    
%If two arguments are passed    
elseif nargin==2
    pZ=varargin{1};
    dimensions=varargin{2};
    
%If too many arguments are passed
else
    disp('Error: too many inputs');
    pZ=varargin{1};
    dimensions=varargin{2};    
end
    
%Check if probabilistic zonotope is too big, such that the computation of
%the plot would be too time consuming
[dim,nrOfGen]=size(pZ.pZ.g);

if nrOfGen<=10
    %get distribution along generators
    [x,prob,probUnscaled]=genDiracs(pZ.pZ,5,1e5);
    %translate points along generators resulting in points and distribution
    %of the xy-plane
    [xy,probXY]=movePoints(pZ.pZ,x,prob);
    %[xy,probXYunscaled]=movePoints(pZ.pZ,x,probUnscaled);
    %plot each point to check correctness of mesh plot
    figure
    colormap('jet')
    
%     plot3(xy(1,:),xy(2,:),probXYunscaled,'LineWidth',2);
%     
%     daspect([5 5 1])
%     %axis([-2 1.8 -2 1 0 1.3])
%     xlim([-2 1.8])
%     ylim([-2 1])
%     zlim([0 1.3])    
    
 
    [x,y,meshProb]=meshNodeProbabilities(pZ,3,30,xy,probXY);
    figure
    mesh(x,y,meshProb);
end
    


%--------------------------------------------------------------------------
function [x,prob,probUnscaled]=genDiracs(pZ,n,maxPoints)
%get dim, Number of generators
[dim,nrOfGen]=size(pZ.g);
%compute infimum and supremum of the generators for the n-sigma bound
var=n*pZ.variance;
infimum=min(pZ.mean-var,[],1);
supremum=max(pZ.mean+var,[],1);    
%for each generator, do the following:
for i=1:nrOfGen
    %define points along generator
    x(i,:) = linspace(infimum(i),supremum(i),maxPoints^(1/nrOfGen));
    y=[];
    %determine prob values for each gaussian distribution
    for j=1:max(find(pZ.variance(:,i))) %do not account for gaussian distribution with zero variance
        y(j,:) = pZ.omega(j,i)*normpdf(x(i,:),pZ.mean(j,i),pZ.variance(j,i));
    end
    %sum up distributions to gaussian mixture distribution
    ySum(i,:)=sum(y,1);
    prob(i,:)=1/sum(ySum(i,:))*ySum(i,:);
    probUnscaled(i,:)=ySum(i,:);
end


%--------------------------------------------------------------------------
function [xy,probXY]=movePoints(pZ,x,prob)
%get dim, Number of generators
[dim,nrOfGen]=size(pZ.g);
%Initialize with points and probabilities of the first generator
xy=pZ.g(:,1)*x(1,:);
probXY=prob(1,:);
%move points along next generator
for i=2:nrOfGen
    %determine new points/probabilities
    tempXY=[];
    tempProbXY=[];
    for j=1:length(prob(1,:))
        %translation matrix TM
        TM=x(i,j)*pZ.g(:,i)*ones(1,length(xy(1,:)));
        tempXY(:,(end+1):(end+length(xy(1,:))))=xy+TM;

        %probability vector
        tempProbXY(:,(end+1):(end+length(xy(1,:))))=prob(i,j)*probXY;
    end
    xy=tempXY;
    probXY=tempProbXY;
end


%--------------------------------------------------------------------------
function [x,y,meshProb]=meshNodeProbabilities(pZ,n,partitions,xy,probXY)

%determine mesh size by n-sigma hyperbox
Z=nSigma(pZ,n);
IH=interval(Z);
inf=infimum(IH);
sup=supremum(IH);
x=linspace(inf(1),sup(1),partitions);
y=linspace(inf(2),sup(2),partitions);
%determine mesh distances
deltaX=x(2)-x(1);
deltaY=y(2)-y(1);
%determine cell volume
A=deltaX*deltaY;
%check total points/ total probability
totalPoints=0;
totalProb=0;

%Build mesh probability matrix
for i=1:length(x)
    %get points that are in the corrext x-bound 
    leftBoundX=x(i)-0.5*deltaX;
    rightBoundX=x(i)+0.5*deltaX;
    inBoundX=(leftBoundX<=xy(1,:)) & (xy(1,:)<rightBoundX);
    xIndices=find(inBoundX);
    %from thesepoints, check for points in the correct y-bound
    for j=1:length(y)
        leftBoundY=y(j)-0.5*deltaY;
        rightBoundY=y(j)+0.5*deltaY;      
        inBoundY=(leftBoundY<=xy(2,xIndices)) & (xy(2,xIndices)<rightBoundY);
        yIndices=find(inBoundY);
        xyIndices=xIndices(yIndices);
        if isempty(xyIndices)
            meshProb(j,i)=0;
            disp('Error! No point in a cell');
        else
            meshProb(j,i)=sum(probXY(xyIndices))/A;
            if meshProb(j,i)>0.1
                disp('problem');
            end
            totalPoints=totalPoints+length(xyIndices);
            totalProb=totalProb+sum(probXY(xyIndices));
        end
    end
end

totalPoints
totalProb
totalProb2=sum(sum(meshProb))*A;
%------------- END OF CODE --------------