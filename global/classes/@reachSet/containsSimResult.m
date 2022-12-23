function res = containsSimResult(R,simRes,varargin)
% containsSimResult - checks if reachable set contains all simulated
%    points of a set of system trajectory; not for hybrid systems
%
% Syntax:  
%    res = containsSimResult(R,simRes)
%    res = containsSimResult(R,simRes,verbose)
%
% Inputs:
%    R - object of class reachSet
%    simRes - object of class simRes
%    verbose - true/false: outputs current simRes.x and total number of
%              simRes.x
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% set default values
verbose = setDefaultValues({false},varargin);

% check input arguments
inputArgsCheck({{R,'att',{'reachSet'},{''}};
                {simRes,'att',{'simResult'},{''}};
                {verbose,'att',{'logical'},{''}}});

% currently not checked cases
if length(R) > 1
    throw(CORAerror('CORA:notSupported','Multiple branches not supported.'));
end


% loop over simulations
for iSim = 1:length(simRes.x)
    
    % output information about status
    if verbose
        disp(['Trajectory ' num2str(iSim) '/' num2str(length(simRes.x)) '...']);
    end
    
    % go over time steps in simulations (skip first point, in initset)
    for timeStep = 2:length(simRes.t{iSim})
        
        % point in trajectory
        simPoint = simRes.x{iSim}(timeStep,:)';
        % time of point
        timePoint = simRes.t{iSim}(timeStep);
           
        % find reachable set at same time as sim point
        for iSet = 1:length(R.timeInterval.time)
            if contains(R.timeInterval.time{iSet},timePoint)
                break;
            end
        end
        
        % check if point in reachable set
        if ~contains(R.timeInterval.set{iSet},simPoint)
            res = false; return;
        end
        
    end
end


% check successful
res = true;

%------------- END OF CODE --------------