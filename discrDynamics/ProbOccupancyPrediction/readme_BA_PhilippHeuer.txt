Changes done by Philipp Heuer in the course of his bachelor thesis:

CORA:

The car class represents a car and consists of the car's segments, the longitudinal position probability
distribution (posProb), the lateral position probability distribution (devProb), the with of the lane of
the car and the segmentLength of the car. Attention: The widthOfLane property does not denote the width
of the lane of the road, but rather the width of the lane of the car. Multiple cars on the same road
can have different widthOfLane-properties. 
Example: car with 5 deviation segments, each 0.5m => widthOfLane = 2.5m

- The conversion from a motion primitive to segments is done in getSegmentsForManeuver.m

- The High-Level Planner introduced in the bachelor thesis can be found in the 'motionplanning' 
  git repository
  
- The Heuristic of the ARA* search is also situated in the 'motionplanning' repository

Old functions that were adapted to be used without the CORA road class:

Old function                    || new function 
intersection.m                  || calculateSingleCrashProbability in calculateCrashProbability.m 
                                || (The function intersection.m has been adapted for an easier use 
                                || with the introduced car class: To calculate the crash probability
                                || between two cars, simply use 
                                || calculateCrashProbability(car1, car2, fArray, 'car');
segPolytope.m                   || getPolytopeForSegment.m 
createPath.m                    || getSegmentsForPath.m
plot.m (in @road)               || plotPositionDistribution.m
segCenter.m                     || segmentCenter.m
segData.m                       || segmentData.m
                                
                                
Other functions:

getMarkovChainResults.m     This method has been created to encapsulate the Markov chain process and 
                            to allow an easy use of the results
                            
How to:
If you want to create a new scenario using the search and the prediction of traffic participants,
please orientate yourself on the scenarios in the 'motionplanning' git repository.                        

Please note: To execute a scenario, CORA as well as the 'motionplanning' git repository have to be checked
out and added to the MATLAB path. Both modules are needed to run the scenarios in the 'motionplanning'
repository (dynamic1simulationPhilipp.m, static1simulationPhilipp.m and static2simulationPhilipp.m).