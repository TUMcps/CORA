function plot3Lanes(obj,except)
% Purpose:  plot road markings of AT08 evading maneuvers
% Pre:      road object
% Post:     ---
% Built:    24.04.08,MA
% Updated:  12.08.09,MA


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

        xLeft(i)=x1-0.5*transLat(1);
        yLeft(i)=y1-0.5*transLat(2);

        xMid1(i)=x1+0.5*transLat(1);
        yMid1(i)=y1+0.5*transLat(2);    
        
        xMid2(i)=x1+1.5*transLat(1);
        yMid2(i)=y1+1.5*transLat(2);           

        xRight(i)=x1+2.5*transLat(1);
        yRight(i)=y1+2.5*transLat(2);     
end

plot(xLeft,yLeft,'k-');
plot(xMid1,yMid1,'k--');
plot(xMid2,yMid2,'k--');
plot(xRight,yRight,'k-');



    
