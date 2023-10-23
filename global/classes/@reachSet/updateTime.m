function R = updateTime(R,time)
% updateTime - updates the times of the reachable set due to the uncertain
%    initial time (currently only supported for one branch in reachSet)
%
% Syntax:
%    R = updateTime(R,time)
%
% Inputs:
%    R - reachSet object
%    time - interval object
%
% Outputs:
%    R - updated reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Authors:       Niklas Kochdumper
% Written:       08-June-2020
% Last update:   18-June-2022 (MW, move from location/reach, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% time shift
deltaT = time - infimum(time);

% shift time of time-point reachable sets
if ~isempty(R.timePoint)
    for i=1:length(R.timePoint.time)
	    R.timePoint.time{i} = R.timePoint.time{i} + deltaT; 
    end
end
% shift time of time-interval reachable sets
if ~isempty(R.timeInterval)
    for i=1:length(R.timeInterval.time)
        R.timeInterval.time{i} = R.timeInterval.time{i} + deltaT;
    end
end

% ------------------------------ END OF CODE ------------------------------
