%converts given maneuver (motion primitive)to segments 
function [segments] = getSegmentsForManeuver(maneuver, segmentLength)
    
    %assemble all x and y values of the current path
    pathX = maneuver.data.x;
    pathY = maneuver.data.y;
    
    fineX = pathX(1):0.1:pathX(end);
    fineY = interp1(pathX, pathY, fineX);   %linearly interpolate the maneuver in small steps of 0.1m
    
    %find points with a distance of approximately segmentLength
    x(1) = fineX(1);    %initial values
    y(1) = fineY(1);
    
    %find interpolated point that has an Euclidean distance to the previous point
    nextIndex = getNextIndexWithApproximateDistance(fineX, fineY, 1, segmentLength);
    while(nextIndex ~= 0)
        x(end+1) = fineX(nextIndex);
        y(end+1) = fineY(nextIndex);
        nextIndex = getNextIndexWithApproximateDistance(fineX, fineY, nextIndex, segmentLength);
    end
    
    %extrapolate the last segment point
    [x, y] = extrapolateEndPoint(x, y, maneuver.data.x(end), maneuver.data.y(end));
    
    angle = calculateAngles(y);
    
    segments.x = x;
    segments.y = y;
    segments.angle = angle;
end

%returns the index of the next point that has a distance that is closest to
%the given distance starting from the point at the given startIndex
function [index] = getNextIndexWithApproximateDistance(x, y, startIndex, distance)
    previousError = inf;
    len = length(x);
    i = startIndex + 1;
    
    while true
        currentDistance = norm([x(startIndex); y(startIndex)] - [x(i); y(i)]);
        currentError = abs(currentDistance - distance);
        
        if(currentError > previousError)
            break
        end
        previousError = currentError;
        i = i + 1;
        
        if(i > len)
            index = 0;
            return
        end
    end
    index = i - 1;
end

%extrapolate endpoint with distance 1 to previous point as a point that
%lies on the straight line spanned by the last point of the given x and y
%vectors and the end points of the maneuver
function [x, y] = extrapolateEndPoint(x, y, maneuverEndX, maneuverEndY)
    dirVector = [maneuverEndX - x(end); maneuverEndY - y(end)];
    %norm directional vector such that the x value is 1
    dirVector = dirVector * (1 / dirVector(1));
    
    %find point on line with approximate distance of 1
    startPoint = [x(end); y(end)];
    factor = 1;
    delta = 0.5;
    while(1)
        potEndPoint = startPoint + (factor * dirVector);
        
        distance = norm([startPoint(1); startPoint(2)] - [potEndPoint(1); potEndPoint(2)]);
        if(abs(distance - 1) > 0.005)
            if(distance > 1)
                factor = factor - delta;
            else
                factor = factor + delta;
            end
            delta = delta / 2;
        else
            break;
        end
    end
    x(end+1) = potEndPoint(1);
    y(end+1) = potEndPoint(2);
end


%caluclate the angles for the given y vector as the y vector determines the angles
function [angles] = calculateAngles(y)
     angles = zeros(1, length(y));

    for i = 1:length(y)-1
        deltaY = y(i+1) - y(i);
        angles(i) = asin(deltaY);
    end
    angles(end) = 0;
end
