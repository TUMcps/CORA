function [tranFrac] = transitionFraction(R)
% transitionFraction - Computes by which probability a transition has
% occured from deceleration to standstill or from acceleration to speed
% limit. The probability is computed straightforward for the time point
% solution; the time interval solution is computed in an approximate manner
% (see AT paper draft)
%
% Syntax:  
%    [tranFrac] = transitionFraction(R,Rcont,tmin,tmax,tend)
%
% Inputs:
%    R - reachable set
%    TP - time struct for initial, last time point in invariant as well as 
%    first and last time point of guard intersection
%
% Outputs:
%    tranFrac - struct for transition probabilities
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      21-April-2009
% Last update:  31-July-2017
%               24-July-2020
% Last revision: ---

%------------- BEGIN CODE --------------

%check if a transition has occured
if length(R)==1
    tranFrac.TP{1} = 1;
    tranFrac.TI{1} = 1;    
else

    % determine the transition fraction for the time point solution;
    % we assume that at maximum only one transition can occur

    % compute volume of reachable set in first location
    vol_1 = volume(R(1).timePoint.set(end));
    %compute volume of reachable set in second location
    vol_2 = volume(R(2).timePoint.set(end));   %use polytope conversion as this has been applied to R, too.
    %compute ratio
    ratio = vol_1/(vol_1 + vol_2);
    
    %set transition fraction for time points
    tranFrac.TP{1} = ratio;
    tranFrac.TP{2} = 1-ratio;
    
    %approximate time interval ratio
    tMax_1 = supremum(R{1}.timePoint.time{end});
    tMax_2 = supremum(R{2}.timePoint.time{end});
    ratioTI = tMax_1/tMax_2;
    
    %set transition fraction for time intervals
    tranFrac.TI{1} = ratioTI;
    tranFrac.TI{2} = 1-ratioTI;    
end


%------------- END OF CODE --------------