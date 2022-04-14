function R = append(R,Radd)
% append - appends one reachSet object at the end of another one
%    (currently, only one branch supported)
%
% Syntax:  
%    obj = add(R,Radd)
%
% Inputs:
%    R - reachSet object
%    Radd - reachSet object
%
% Outputs:
%    R - resulting reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Author:       Mark Wetzlinger
% Written:      18-May-2021         
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% currently, only one branch per reachSet object
if length(R) > 1 || length(Radd) > 1
    error('Multiple branches of R not checked.');
end

% empty objects
if isempty(R)
    R = Radd;
elseif isempty(Radd)
    % just R
else        
    % general case:

    % get final time of R
    shift = R.timePoint.time{end};
    
    % number of sets of Radd
    noSets = length(Radd.timePoint.set);
    
    % time-point sets
    R.timePoint.set = [R.timePoint.set; Radd.timePoint.set];
    % shift time
    for i=1:noSets
        Radd.timePoint.time{i} = Radd.timePoint.time{i} + shift;
    end
    R.timePoint.time = [R.timePoint.time; Radd.timePoint.time];
    
    
    % time-interval sets
    if ~isempty(R.timeInterval.set)
        
        R.timeInterval.set = [R.timeInterval.set; Radd.timeInterval.set];
        % shift time
        for i=1:noSets
            Radd.timeInterval.time{i} = Radd.timeInterval.time{i} + shift;
        end
        R.timeInterval.time = [R.timeInterval.time; Radd.timeInterval.time];
        
    end
    
end

%------------- END OF CODE --------------