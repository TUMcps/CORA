function han = plotCircle(center,radius,NOP,style)
% plotCircle - plots a circle
%
% Syntax:
%    han = plotCircle(center,radius,NOP,style)
%
% Inputs:
%    center - center of circle
%    radius - radius of circle
%    NOP - number of points on the circle
%    style - LineSpec or Name-Value pairs (same as for plot)
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%     plotCircle([1,3],3,1000,':'); 
%     plotCircle([2,4],2,1000,'--');

% Authors:       Zhenhai Wang
% Written:       December-2002
% Last update:   17-June-2022 (MW, only formatting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs',3));
elseif nargin==3
    style = 'b-';
end

THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);

%H=plot(X,Y,style); %<-- modified
if strcmp(style,'b-')
    %H=fill(X,Y,[.8 .8 .8],'EdgeColor','k'); %<-- modified
    han=plot(X,Y,'k'); %<-- modified
else
    han=fill(X,Y,[.95 .95 .95],'EdgeColor','k'); %<-- modified
end

%axis square;

% ------------------------------ END OF CODE ------------------------------
