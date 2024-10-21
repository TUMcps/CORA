classdef tentativeKleeneSignal
% tentativeKleeneSignal - a Kleene signal allowing for tentatively true and tentatively false values
%
% When observing a new value, the new value will be merged with the old value.
% If contradicting values are observed, the signal will change its value to unknown.
%
% Syntax:
%    sig = tentativeKleeneSignal(canBeTrue,canBeFalse)
%
% Inputs:
%    canBeTrue - Boolean pointSegmentSignal indicating whether the value could be true
%    canBeFalse - Boolean pointSegmentSignal indicating whether the value could be false
%
% Outputs:
%    sig - generated tentativeKleeneSignal object
%
% Example:
%    canBeTrue = pointSegmentSignal.indicator(stlInterval(1,3),true,false);
%    canBeFalse = pointSegmentSignal.indicator(stlInterval(2,4),true,false);
%    sig = tentativeKleeneSignal(canBeTrue,canBeFalse);
%    sig.plot()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       14-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    canBeTrue pointSegmentSignal
    canBeFalse pointSegmentSignal
end

methods
    % constructor
    function obj = tentativeKleeneSignal(canBeTrue,canBeFalse)
        obj.canBeTrue = canBeTrue;
        obj.canBeFalse = canBeFalse;
    end

    % set a new tentative value in the given interval
    % if the new value contradicts the old value, the signal will be set to unknown
    function sig = observe(obj,interval,value)
        switch value
            case kleene.True
                canBeTrueSig = obj.canBeTrue.set(interval,true);
                canBeFalseSig = obj.canBeFalse;
                case kleene.False
                canBeTrueSig = obj.canBeTrue;
                canBeFalseSig = obj.canBeFalse.set(interval,true);
            case kleene.Unknown
                canBeTrueSig = obj.canBeTrue.set(interval,true);
                canBeFalseSig = obj.canBeFalse.set(interval,true);
            otherwise
                assert(false,'Not a known Kleene value');
        end

        sig = tentativeKleeneSignal(canBeTrueSig,canBeFalseSig);
    end

    % remove all observations inside the interval
    function sig = setInvalid(obj,interval)
        canBeTrueSig = obj.canBeTrue.set(interval,false);
        canBeFalseSig = obj.canBeFalse.set(interval,false);
        sig = tentativeKleeneSignal(canBeTrueSig,canBeFalseSig);
    end

    % extract the observations inside the interval
    function sig = getInterval(obj,interval)
        if isemptyobject(interval)
            % everything is outside the interval and thus invalid
            sig = fourValuedSignal.uniformSignal(fourValued.Inconclusive);
            return;
        end

        % invalidate outside the interval
        left = interval.toLeft();
        right = interval.toRight();
        obj = obj.setInvalid(left);
        obj = obj.setInvalid(right);
        
        sig = obj.toFourValued();
    end

    % plot the signal
    function han = plot(obj)
        if nargout > 0
            han = plot(obj.toFourValued());
        else
            plot(obj.toFourValued());
        end
    end  
end

methods (Access = private)
    % convert the current observations to a four-valued signal
    function sig = toFourValued(obj)
        unk = obj.canBeTrue & obj.canBeFalse;
        kleeneSig = kleeneSignal(obj.canBeTrue,unk);
        inc = ~(obj.canBeTrue | obj.canBeFalse);
        sig = fourValuedSignal(kleeneSig,inc);
    end
end

methods (Static)
    % create an empty signal that has not yet observed any values
    function obj = emptySignal
        canBeTrueSig = pointSegmentSignal.indicator(stlInterval(0,inf),false,false);
        canBeFalseSig = pointSegmentSignal.indicator(stlInterval(0,inf),false,false);
        obj = tentativeKleeneSignal(canBeTrueSig,canBeFalseSig);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
