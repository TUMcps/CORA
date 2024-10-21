classdef pointSegmentSignal < logicSignal
% pointSegmentSignal - a continuous-time signal allowing values to hold at singular points in time
%
% The timePoints array must contain strictly increasing values and start at 0.
% The values array must contain a value for each time point and its following interval,
% i.e. it has twice the length of timePoints.
% No time point may have the same value as both the preceding and following intervals.
%
% Syntax:
%    sig = pointSegmentSignal(timePoints,values)
%
% Inputs:
%    timePoints - strictly increasing array of time points
%    values - array of values, twice the length of timePoints (odd indices are values at time points, even indices are values at intervals)
%
% Outputs:
%    sig - generated pointSegmentSignal object
%
% Example:
%    timePoints = [0,1,2];
%    values = [true,false,false,true,true,false];
%    sig = pointSegmentSignal(timePoints,values);
%    sig.plot()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: logicSignal

% Authors:       Florian Lercher
% Written:       12-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    timePoints
    values
end

methods
    % constructor
    function sig = pointSegmentSignal(timePoints,values)
        if CHECKS_ENABLED
            aux_checkInputs(timePoints,values);
        end
        sig.timePoints = timePoints;
        sig.values = values;
    end

    % methods in other files
    res = allTrue(obj,interval,cond)
    res = anyTrue(obj,interval,cond)
    sig = cutAtFirstFallingEdge(obj)
    intervals = findIntervals(obj,cond)
    [time,inclusive] = findLast(obj,value)
    han = plot(sig)
    sig = set(obj,interval,value)
    sig = until(lhs,interval,rhs)

    % signal value at a given time
    function val = at(obj,time)
        val = priv_atIdx(obj,time);
    end

    % logical negation
    function sig = not(obj)
        sig = pointSegmentSignal(obj.timePoints,~obj.values);
    end
    
    % signal equality
    function equal = eq(obj,other)
        if isa(other,'pointSegmentSignal')
            equal = all(withinTol(obj.timePoints,other.timePoints,eps)) && isequal(obj.values,other.values);
        else
            equal = false;
        end
    end
end

methods (Static)
    % methods in other files
    sigs = combine(op,varargin)
    sig = indicator(interval,val,default)

    % logical conjunction
    function sig = and_(varargin)
        sig = pointSegmentSignal.combine(@all,varargin{:});
    end

    % logical disjunction
    function sig = or_(varargin)
        sig = pointSegmentSignal.combine(@any,varargin{:});
    end
end

methods (Access = private)
    % discrete length of the signal
    function l = length(obj)
        l = length(obj.timePoints);
    end

    % point value for a given index into the timePoints array
    function val = pointValue(obj,i)
        val = obj.values(2 * i - 1);
    end

    % preceding interval value for a given index into the timePoints array
    function val = precIntervalValue(obj,i)
        val = obj.values(2 * i - 2);
    end

    % succeeding interval value for a given index into the timePoints array
    function val = succIntervalValue(obj,i)
        val = obj.values(2 * i);
    end

    % preceding time point for a given index into the values array
    function [time,atPoint] = precTimeAtValIdx(obj,i)
        ti = floor((i + 1) / 2);
        time = obj.timePoints(ti);
        atPoint = mod(i,2) == 1;
    end

    % succeeding time point for a given index into the values array
    function [time,atPoint] = succTimeAtValIdx(obj,i)
        ti = floor(i / 2) + 1;
        if ti > length(obj.timePoints)
            time = inf;
        else
            time = obj.timePoints(ti);
        end
        atPoint = mod(i,2) == 1;
    end
end
end


% Auxiliary functions -----------------------------------------------------

function aux_checkInputs(timePoints,values)
    % signal may not be empty
    if isempty(timePoints) || isempty(values)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'The time and value arrays must not be empty.'));
    end

    % for each time point there must be a value for the point itself and the following interval
    if 2 * length(timePoints) ~= length(values)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'There must be a value for each time point and its following interval.'));
    end

    % the initial time point must be 0
    if timePoints(1) ~= 0
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'The first time point must be 0.'));
    end

    % points in time must be strictly increasing
    if ~issorted(timePoints, 'strictascend')
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'The time array must be sorted.'));
    end

    % the value at a time point may not match the values of both the preceding and the following interval
    for i=2:length(timePoints)
        valI = 2 * i - 1;
        if values(valI) == values(valI + 1) && values(valI) == values(valI - 1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'The value at a time point may not match both the preceding and following interval.'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
