%calculates the crashProbability between the two given cars
function [crashProb] = calculateCrashProbability(car1, car2, fArray, type)
    %tic
    segments1 = get(car1, 'segments');
    posProb1 = get(car1, 'posProb');
    devProb1 = get(car1, 'devProb');
    widthOfLane1 = get(car1, 'widthOfLane');

    segments2 = get(car2, 'segments');
    posProb2 = get(car2, 'posProb');
    devProb2 = get(car2, 'devProb');
    widthOfLane2 = get(car2, 'widthOfLane');

    %Check inputs
    if(length(posProb1) ~= length(posProb2))
        assert(length(posProb1) == 1 || length(posProb2) == 1);
    end
    
    %in case of one car being static with only one scalar as posProb
    maxLengthPosProb = max(length(posProb1), length(posProb2));
    
    %make sure both posProb have the same length
    posProb1 = convertToCellArrayWithLength(posProb1, maxLengthPosProb);
    posProb2 = convertToCellArrayWithLength(posProb2, maxLengthPosProb);

    %space pre-allocation
    crashProb = zeros(1, maxLengthPosProb);

    for iStep = 1:maxLengthPosProb
        crashProb(iStep) = calculateSingleCrashProbability(segments1, segments2, posProb1{iStep}, posProb2{iStep}, ...
            devProb1, devProb2, widthOfLane1, widthOfLane2, fArray, type);
    end
    %toc
end

%returns the given posProb as a cell array with the given length since 
%posProbs can be either already a cell array, a single cell or a scalar
function cells = convertToCellArrayWithLength(posProb, neededLength) 
    if(length(posProb) == 1)
        cells = cell(1, neededLength);
        if(iscell(posProb))
            cells(:) = posProb;
        else
            cells(:) = {posProb};
        end
    else
        cells = posProb;
    end
end

%Adapted function intersection without use of CORA road class
function [p] = calculateSingleCrashProbability(segments1, segments2, posProb1, posProb2, devProb1, devProb2, widthOfLane1, widthOfLane2, fArray, type)
    % intersection - computes the probability that two reachable sets 
    % of traffic participants intersect.
    %
    % Inputs:
    %    segments1 - path segments of traffic participant 1
    %    segments2 - path segments of traffic participant 2
    %    posProb1 - segment probability vector of TP 1
    %    posProb2 - segment probability vector of TP 2
    %    devProb1 - deviation probability distribution of TP 1
    %    devProb2 - deviation probability distribution of TP 2
    %    widthOfLane1 - width of the lane of TP1
    %    widthOfLane2 - width of the lane of TP2
    %    fArray - result of intersectionDatabase.m
    %    type - type of the vehicle (options are 'car', 'carCons' and 'bicycle')
    %
    % Outputs:
    %    p - intersection probability

    %------------- BEGIN CODE --------------

    if strcmp(type,'car') || strcmp(type,'carCons')
        %obtain radii of enclosing circles of road segments
        r = fArray.radius.r.car;
        rComb = fArray.radius.rComb.car;
        %load intersection database
        intDatabase = fArray.val.car;
        if strcmp(type, 'carCons') % choose between normal and conservative intersection
            intFlag = 1;
        else 
            intFlag = 0;
        end
    elseif strcmp(type,'bicycle')
        %obtain radii of enclosing circles of road segments
        r = 0.5 * (fArray.radius.r.car + fArray.radius.r.bicycle);
        rComb = 0.5 * (fArray.radius.rComb.car + fArray.radius.rComb.bicycle);
        %load intersection database
        intDatabase = fArray.val.bicycle;    
        if strcmp(type,'carCons') % choose between normal and conservative intersection
            intFlag = 1;
        else 
            intFlag = 0;
        end
    end    

    %compute segment center vectors of R1
    ind1 = find(posProb1 > 0);
    [c1] = cellCenter(segments1, ind1, [], widthOfLane1, length(devProb1));

    %compute segment center vectors of R2
    ind2 = find(posProb2 > 0);
    [c2] = cellCenter(segments2, ind2, [], widthOfLane2, length(devProb2));

    %init potPairs
    potPairs = [];

    %obtain potential segment intersection pairs
    if (~isempty(c1) & ~isempty(c2))
        for i = 1:length(c1(:,1))
            for j = 1:length(c2(:,1))
                dist = norm(c1(i,:) - c2(j,:));
                if dist < rComb
                    potPairs(end+1,:) = [ind1(i), ind2(j)];
                end
            end
        end
    end

    %init probability
    p = 0;

    if ~isempty(potPairs)
        for iPair = 1:length(potPairs(:,1))
            %compute deviation numbers whose polytopes might intersect
            [potDevPairs] = devNumbers(segments1, segments2, potPairs(iPair,:), r, devProb1, devProb2, widthOfLane1, widthOfLane2);

            if ~isempty(potDevPairs)
                for k = 1:length(potDevPairs(:,1))
                    %generate polytopes
                    [pos1, angle1] = segmentData(segments1, potPairs(iPair,1), potDevPairs(k,1), length(devProb1), widthOfLane1);
                    [pos2, angle2] = segmentData(segments2, potPairs(iPair,2), potDevPairs(k,2), length(devProb2), widthOfLane2);

                    intersected = intersection_db(fArray, intDatabase, pos1, angle1, pos2, angle2);
                    
                    if intersected > 0
                        %compute partial probabilities
                        prob1 = posProb1(potPairs(iPair,1)) * devProb1(potDevPairs(k,1));
                        prob2 = posProb2(potPairs(iPair,2)) * devProb2(potDevPairs(k,2));

                        %add to probability
                        if intFlag
                            p = p + (prob1 * prob2);
                        else
                            p = p + (intersected * prob1 * prob2);
                        end
                    end
                end
            end
        end
    end
end


%auxiliary function to get pairs of segments that potentially intersect
function [potPairs] = devNumbers(segments1, segments2, potPair, r, devProb1, devProb2, widthOfLane1, widthOfLane2)

    %get deviation segments with nonzero probability
    ind1 = find(devProb1);
    ind2 = find(devProb2);

    %compute centers of deviation polytopes
    [cDev1] = cellCenter(segments1, potPair(1), ind1, widthOfLane1, length(devProb1));
    [cDev2] = cellCenter(segments2, potPair(2), ind2, widthOfLane2, length(devProb2));

    %init potPairs
    potPairs = [];

    %obtain potential deviation polytope intersection pairs
    for i = 1:length(cDev1(:,1))
        for j = 1:length(cDev2(:,1))
            dist = norm(cDev1(i,:) - cDev2(j,:));
            if dist < r
                potPairs(end+1,:) = [ind1(i), ind2(j)];
            end
        end
    end
end

%adapted function intersection_database without use of CORA road class
function p = intersection_db(fArray,intDatabase,pos1,angle1,pos2,angle2)
% intersection_db - checks whether two rectangles intersect based on
% a database; an intersection probability is provided for a given
% uncertainty of the center of the rectangles; The database lookup results
% in a small quantization error

    %nr of segments
    nrOfAngleSegments=length(intDatabase);
    [nrOfxSegments,nrOfySegments]=size(intDatabase{1});

    %generate rotation matrix
    rot=[cos(-angle1) -sin(-angle1);...
         sin(-angle1) cos(-angle1)];
    %return relAngle
    relAngle=angle2-angle1;
    relAngle=rem(relAngle,pi);
    %return relPosition
    relPos=rot*(pos2-pos1);
    %generate positive position differences
    if prod(relPos)<0
        relAngle=-relAngle;
    end
    relPos=abs(relPos);
    %generate positive angle differences    
    if relAngle<0
        relAngle=relAngle+pi;
    end   

    %obtain angle and position segments
    iAngleSeg=round(relAngle/fArray.segLength.angle)+1;
    iXseg=round(relPos(1)/fArray.segLength.x)+1;
    iYseg=round(relPos(2)/fArray.segLength.y)+1;

    %in case some segments are zero
    if (iAngleSeg==0) || (iAngleSeg>nrOfAngleSegments)
        p = 0;
        return;
    end
    if (iXseg==0) || (iXseg>nrOfxSegments)
        p = 0;
        return;
    end   
    if (iYseg==0) || (iYseg>nrOfySegments) 
        p = 0;
        return;
    end  

    %look up intersection probability
    p = intDatabase{iAngleSeg}(iXseg,iYseg);
end
