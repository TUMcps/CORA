classdef kleeneSignalBuilder < matlab.mixin.Copyable
% kleeneSignalBuilder - construct a three-valued signal of Kleene logic by
% defining the time intervals where a property can be true or false
%
% Syntax:
%    builder = kleeneSignalBuilder(10)
%
% Inputs:
%    dur - duration of the signal
%
% Outputs:
%    obj - kleeneSignalBuilder object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: kleene signal

% Authors:       Benedikt Seidl
% Written:       23-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % boolean signal where the property can be true
    possiblyTrueSignal signal

    % boolean signal where the property can be false
    possiblyFalseSignal signal
end

methods

    function obj = kleeneSignalBuilder(dur)
        obj.possiblyTrueSignal = signal(dur, false);
        obj.possiblyFalseSignal = signal(dur, false);
    end

    function setPossiblyTrue(obj, int)
        % obtain duration of signal
        dur = duration(obj.possiblyTrueSignal);

        % construct indicator signal for given interval
        ind = signal.indicator(dur, int, true);

        % combine signals
        obj.possiblyTrueSignal = obj.possiblyTrueSignal | ind;
    end

    function setPossiblyFalse(obj, int)
        % obtain duration of signal
        dur = duration(obj.possiblyFalseSignal);

        % construct indicator signal for given interval
        ind = signal.indicator(dur, int, true);

        % combine signals
        obj.possiblyFalseSignal = obj.possiblyFalseSignal | ind;
    end

    function sig = kleeneSignal(obj)
        % construct three-valued signal from boolean indicator signals
        sig = combine(obj.possiblyTrueSignal, ...
            obj.possiblyFalseSignal, @comb);

        function out = comb(possiblyTrue, possiblyFalse)
            if possiblyTrue && ~possiblyFalse
                % property can only be true
                out = kleene.True;
            elseif possiblyFalse && ~possiblyTrue
                % property can only be false
                out = kleene.False;
            elseif possiblyTrue && possiblyFalse
                % property can be both true and false
                out = kleene.Unknown;
            else
                % With masks a property can be neither true nor false. We
                % set the signal to Unknown for convenience.
                out = kleene.Unknown;
            end
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
