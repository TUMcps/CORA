classdef finiteSignal < logicSignal
% finiteSignal - a continuous-time discrete-valued signal without Zeno behaviour
%
% The time and value arrays must have the same length. The time array must
% be sorted in ascending order and no singleton intervals are allowed,
% i.e., no value may hold only at a singular point in time. Also, two
% consecutive entries in the value array may contain the same value.
%
% Syntax:
%    sig = finiteSignal(time, value)
%
% Inputs:
%    time - points in time until the corresponding values are valid
%    value - value of the signal at the corresponding interval
%
% Outputs:
%    sig - generated signal object
%
% Example:
%    tt = kleene.True; ff = kleene.False; uu = kleene.Unknown;
%
%    time = [1.2 2.4 2.7 3.9 4.0 4.5 6.3];
%    value = [tt uu ff uu ff tt uu];
%
%    sig = finiteSignal(time, value);
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       11-August-2022
% Last update:   08-February-2024 (FL, rename from signal to finiteSignal)
%                09-February-2024 (FL, implement logicSignal interface)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % points in time until the corresponding values are valid
    time {mustBeNonnegative}

    % value of the signal at the corresponding interval
    value
end

methods

    % class constructor
    function sig = finiteSignal(time, value)
        % signal may not be empty
        if isempty(time) || isempty(value)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'The time and value arrays must not be empty.'));
        end

        % length of arrays must match
        if length(time) ~= length(value)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'The time and value arrays must have same length.'));
        end

        % points in time must be strictly increasing
        if ~issorted(time, 'strictascend')
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'The time array must be sorted.'));
        end

        % no identical consecutive values
        for i=1:length(value)-1
            if value(i) == value(i+1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'No identical consecutive values allowed.'));
            end
        end

        sig.time = time;
        sig.value = value;
    end

    % negation
    function out = not(sig)
        for i=1:length(sig)
            val(i) = ~ sig.value(i);
        end

        out = finiteSignal(sig.time, val);
    end

    % length
    function out = length(sig)
        out = length(sig.time);
    end

    % duration
    function out = duration(sig)
        out = sig.time(end);
    end

    val = at(sig, time)
    sig = set(obj, interval, value)
    sig = until(lSig, int, rSig)
end

methods (Static)

    % indicator signal for the given interval
    sig = indicator(dur, int, val)

    % point-wise combination
    sigs = combine(op, varargin)

    % conjunction
    function out = and_(varargin)
        out = finiteSignal.combine(@and, varargin{:});
    end

    % disjunction
    function out = or_(varargin)
        out = finiteSignal.combine(@or, varargin{:});
    end

end
end

% ------------------------------ END OF CODE ------------------------------
