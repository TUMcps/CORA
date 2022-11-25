function [segments] = getSegmentsForPath(initials, deltaAngles, timeSteps, segmentLength)
%replaces createPath.m without using the CORA road class
    plotPath = 0;

    %Preallocation of space
    lengthOfArrays = sum(timeSteps) + 1;
    angle = zeros(1, lengthOfArrays);
    x = zeros(1, lengthOfArrays);
    y = zeros(1, lengthOfArrays);

    %initial values for angle, x and y
    angle(1) = initials(1);
    x(1) = initials(2);
    y(1) = initials(3);
    j = 1;

    for i = 1:(length(deltaAngles))
        for iStep = 1:timeSteps(i)
            angle(j+1) = angle(j) + deltaAngles(i);
            %get x-values
            x(j+1) = x(j) + segmentLength * cos(angle(j));
            %get y-values
            y(j+1) = y(j) + segmentLength * sin(angle(j));    
            %counter
            j = j + 1;
        end
    end

    %plot the path
    if plotPath
        figure
        plot(x,y, '-x');
        axis equal
    end

    %write to object
    segments.x = x;
    segments.y = y;
    segments.angle = angle;

    save('path.mat', 'x', 'y', 'angle');
end
