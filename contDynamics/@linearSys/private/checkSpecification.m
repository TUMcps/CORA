function [res,YtimeInt,YtimePoint] = checkSpecification(spec,XtimeInt,YtimeInt,YtimePoint,idx)
% checkSpecification - check safety properties for current time-interval
%    reachable set; if a violation occurs, return truncated structs
%
% Syntax:
%    [res,timeInt,timePoint] = checkSpecification(spec,timeInt,timePoint,idx)
%
% Inputs:
%    spec - object of class specification
%    XtimeInt - struct about time-interval reachable set
%    YtimeInt - struct about time-interval output set
%    YtimePoint - struct about time-point output set
%    idx - index of current time interval 
%
% Outputs:
%    res - true if specifications are satisfied, otherwise false
%    YtimeInt - (truncated) struct about time-interval output set
%               (only required if a violation was detected)
%    YtimePoint - (truncated) struct about time-point output set
%               (only required if a violation was detected)
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-November-2022
% Last update:   07-December-2022 (MW, bug fix for spec.type = 'invariant')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init satisfaction
res = true;

for i=1:length(spec)
    if strcmp(spec(i).type,'invariant')
        % specification is an invariant -> called in hybrid system analysis to
        % check if the reachable set has left the invariant, therefore we use
        % the reachable set instead of the output set
        
        % the linearSys algorithms 'decomp', 'krylov', and 'adaptive' cannot
        % be used in conjunction with hybrid system analysis; thus, the call to
        % this function has XtimeInt = []
        if isnumeric(XtimeInt) && isempty(XtimeInt)
            throw(CORAerror('CORA:notSupported',...
                sprintf(['The chosen algorithm for linear systems (options.linAlg)\n',...
                'cannot be used for the analysis of hybrid systems.'])));
        end
    
        if ~check(spec(i),XtimeInt,YtimeInt.time{idx})
            % violation
            res = false;
            % truncate reachable set until current time interval
            YtimeInt.set = YtimeInt.set(1:idx);
            YtimeInt.time = YtimeInt.time(1:idx);
            % index for time-point shifted by one as initial set at index 1
            YtimePoint.set = YtimePoint.set(1:idx+1);
            YtimePoint.time = YtimePoint.time(1:idx+1);
            return
        end
    else
        % specification on the output set
        if ~check(spec(i),YtimeInt.set{idx},YtimeInt.time{idx})
            % violation
            res = false;
            % truncate output set until current time interval
            YtimeInt.set = YtimeInt.set(1:idx);
            YtimeInt.time = YtimeInt.time(1:idx);
            % index for time-point shifted by one as initial set at index 1
            YtimePoint.set = YtimePoint.set(1:idx+1);
            YtimePoint.time = YtimePoint.time(1:idx+1);
            return
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
