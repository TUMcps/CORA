function example_probOccupancyPrediction_carReachCircles()
% example_probOccupancyPrediction_carReachCircles - ???
%
% Syntax:
%    example_probOccupancyPrediction_carReachCircles()
%
% Inputs:
%    ---       
%
% Outputs:
%    ---
% 
% Example:
%    carReachCircles();

% Authors:       Matthias Althoff
% Written:       19-August-2009
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%unknown plan

%plot road 
R=road(4,5,7); % width:4, seglength:1
R=createPath(R,[0,-10,0],[0],[13]);
plotCrossing(R,[13,13]);
hold on

% figure;
% hold on

for iCar=1:2
%for iCar=1:1

    if iCar==1
        %set parameters
        mu=1; 
        g=9.81; %[m/s^2]
        v0=20; %[m/s]
        %v0=10; %[m/s]
        p0=[0;0]; %[m]
        %p0=[10;0]; %[m]
        plotType='b-';
    else
        %set parameters
        mu=1; 
        g=9.81; %[m/s^2]
        v0=-20; %[m/s]
        p0=[30;4]; %[m]     
        plotType='g-';
    end

    %set time points
    %t=0:0.02:1;
    %t=0:0.1:1;
    t=0:0.02:1;

    %compute radii and centers
    r=0.5*mu*g*t.^2;
    c=v0*t;

    %plot circles
    for i=1:length(t)
        plotCircle(p0+[c(i);0],r(i),100,plotType);
    end
    
%     %plot centers
%     for i=1:length(t)
%         plot(p0(1)+c(i),p0(2),'k+');
%     end
end

axis equal


% %--------------------------------------------------------------------------
% %known plan
% 
% %plot road 
% R=road(4,5,7); % width:4, seglength:1
% R=createPath(R,[0,-10,0],[0],[13]);
% plotCrossing(R,[13,13]);
% hold on
% 
% % figure;
% % hold on
% 
% for iCar=1:2
% %for iCar=1:1
% 
%     if iCar==1
%         %set parameters
%         mu=1; 
%         g=9.81; %[m/s^2]
%         v0=20; %[m/s]
%         %v0=10; %[m/s]
%         p0=[0;0]; %[m]
%         %p0=[10;0]; %[m]
%         %plotType='b-';
%         plotType='k-';
%     else
%         %set parameters
%         mu=1; 
%         g=9.81; %[m/s^2]
%         v0=-20; %[m/s]
%         p0=[40;4]; %[m]     
%         %plotType='g-';
%         plotType='k-';
%     end
% 
%     %set time points
%     %t=0:0.02:1;
%     %t=0:0.1:1;
%     t=0:0.02:1;
% 
%     %compute radii and centers
%     c=v0*t;
% 
%     %plot circles
%     w=2.5;
%     l=5.5;
%     deltaVec=[0.5*l;0.5*w];
%     for i=1:length(t)
%         center=p0+[c(i);0];
%         I=interval([center-deltaVec,center+deltaVec]);
%         plot(I,[1 2],'grayFrame');
%     end
% end
% 
% axis equal

% ------------------------------ END OF CODE ------------------------------
