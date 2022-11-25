function plotCrossing(obj,except)
% Purpose:  plot road markings of the IV08 crossing
% Pre:      road object
% Post:     ---
% Built:    22.11.07,MA


hold on

%load road variables
x=obj.segments.x;
y=obj.segments.y;
angle=obj.segments.angle;

for i=1:except(1)
    
        %obtain x and y coordinates
        x1=x(i);
        y1=y(i);
        angle1=angle(i);

        transLat(1)=cos(angle1-0.5*pi)*1*obj.width;
        transLat(2)=sin(angle1-0.5*pi)*1*obj.width;

        xLeft(i)=x1-1.5*transLat(1);
        yLeft(i)=y1-1.5*transLat(2);

        xMid(i)=x1-0.5*transLat(1);
        yMid(i)=y1-0.5*transLat(2);    

        xRight(i)=x1+0.5*transLat(1);
        yRight(i)=y1+0.5*transLat(2);     
end

plot(xLeft,yLeft,'k-');
plot(xMid,yMid,'k--');
plot(xRight,yRight,'k-');


xLeft=[]; yLeft=[];
xMid=[]; yMid=[];
xRight=[]; yRight=[];

for i=1:length(x)-except(2)
    
        %obtain x and y coordinates
        x1=x(except(2)+i);
        y1=y(except(2)+i);
        angle1=angle(except(2)+i);

        transLat(1)=cos(angle1-0.5*pi)*1*obj.width;
        transLat(2)=sin(angle1-0.5*pi)*1*obj.width;

        xLeft(i)=x1-1.5*transLat(1);
        yLeft(i)=y1-1.5*transLat(2);

        xMid(i)=x1-0.5*transLat(1);
        yMid(i)=y1-0.5*transLat(2);    

        xRight(i)=x1+0.5*transLat(1);
        yRight(i)=y1+0.5*transLat(2);     
end

plot(xLeft,yLeft,'k-');
plot(xMid,yMid,'k--');
plot(xRight,yRight,'k-');
    
