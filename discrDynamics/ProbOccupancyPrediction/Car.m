classdef Car < matlab.mixin.SetGet
    
    properties
        segments
        posProb
        devProb
        widthOfLane %ATTENTION: widthOfLane does not denote the width of the lane of the road, but rather the width of the lane of the car
                    %Multiple cars on the same road can have different widthOfLane depending on the width of lateral uncertainty
                    %Example: car with 5 deviation segments, each 0.5m => widthOfLane = 2.5
        segmentLength;
    end
    
    methods
        function obj = Car(segments, posProb, devProb, widthOfLane, segmentLength)
            obj.segments = segments;
            obj.posProb = posProb;
            obj.devProb = devProb;
            obj.widthOfLane = widthOfLane;
            obj.segmentLength = segmentLength;
        end
    end
    
end

