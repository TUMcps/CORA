function sigs = combine(op,varargin)
% combine - compute the point-wise combination of signals using the given operator
%
% If n = length(varargin), op has to process an 1xn array of values.
% If op returns a 1xm array of values, the output sigs is an 1xm array of signals.
%
% Syntax:
%    sigs = pointSegmentSignal.combine(op,sig1,sig2,...)
%
% Inputs:
%    op - function handle to the operator
%    varargin - pointSegmentSignal objects
%
% Outputs:
%    sigs - array of result signals
%
% Example:
%    sig1 = pointSegmentSignal([0 2 3],[true true false false false true]);
%    sig2 = pointSegmentSignal([0 1 3],[true true false false false true]);
%    sig = pointSegmentSignal.combine(@all,sig1,sig2);
%    assert(sig == sig2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       12-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

numInputs = length(varargin);
numOutputs = length(op(false(1,numInputs)));
inIdx = ones(1,numInputs);
outIdx = ones(1,numOutputs);
tp = cell(1,numOutputs);
val = cell(1,numOutputs);
while any(inIdx <= cellfun(@length,varargin))
    times = arrayfun(@aux_time,varargin,inIdx);
    timePoint = min(times);
    isMin = times - timePoint < eps;
    pointValues = op(arrayfun(@aux_pv,varargin,inIdx,isMin));
    intValues = op(arrayfun(@aux_iv,varargin,inIdx,isMin));
    minIdx = find(isMin);
    inIdx(minIdx) = inIdx(minIdx) + 1;

    for i = 1:numOutputs
        pointValue = pointValues(i);
        intValue = intValues(i);
        j = outIdx(i);
        if j > 1
            % check whether the new point is equal to the preceding and following interval
            if pointValue == val{i}(2 * j - 2) && pointValue == intValue
                % if so, skip this point
                continue;
            end
        end
        tp{i}(j) = timePoint;
        val{i}(2 * j - 1) = pointValue;
        val{i}(2 * j) = intValue;
        outIdx(i) = j + 1;
    end
end

for i = 1:numOutputs
    sigs(i) = pointSegmentSignal(tp{i},val{i});
end
end


% Auxiliary functions -----------------------------------------------------

function val = aux_time(signal,idx)
    if idx <= length(signal{1}.timePoints)
        val = signal{1}.timePoints(idx);
    else
        val = inf;
    end
end

function val = aux_pv(signal,idx,atPoint)
    if atPoint
        val = signal{1}.pointValue(idx);
    else
        val = signal{1}.precIntervalValue(idx);
    end
end

function val = aux_iv(signal,idx,atPoint)
    if atPoint
        val = signal{1}.succIntervalValue(idx);
    else
        val = signal{1}.precIntervalValue(idx);
    end
end

% ------------------------------ END OF CODE ------------------------------
