function R = append(R,Radd)
% append - appends one reachSet object at the end of another one
%    (currently, only one branch supported)
%
% Syntax:
%    obj = append(R,Radd)
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

% Authors:       Mark Wetzlinger
% Written:       18-May-2021         
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% currently, only one branch per reachSet object
if length(R) > 1 || length(Radd) > 1
    throw(CORAerror('CORA:notSupported','Multiple branches not supported.'));
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
    nrSets_tp = length(Radd.timePoint.set);
    nrSets_ti = length(Radd.timeInterval.set);
    
    % time-point sets
    R.timePoint.set = [R.timePoint.set; Radd.timePoint.set];
    % shift time
    for i=1:nrSets_tp
        Radd.timePoint.time{i} = Radd.timePoint.time{i} + shift;
    end
    R.timePoint.time = [R.timePoint.time; Radd.timePoint.time];
    
    
    % time-interval sets
    if ~isempty(R.timeInterval.set)
        
        R.timeInterval.set = [R.timeInterval.set; Radd.timeInterval.set];
        % shift time
        for i=1:nrSets_ti
            Radd.timeInterval.time{i} = Radd.timeInterval.time{i} + shift;
        end
        R.timeInterval.time = [R.timeInterval.time; Radd.timeInterval.time];
        
    end
    
end

% ------------------------------ END OF CODE ------------------------------
