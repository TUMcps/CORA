classdef fourValuedSignal < logicSignal
% fourValuedSignal - a signal over four-valued truth values represented by a Kleene and a Boolean signal
%
% Since the four-valued signal is represented by a fourValuedSignal and a pointSegmentSignal object,
% most operations are simple wrappers.
%
% Syntax:
%    sig = fourValuedSignal(signal,isInconclusive)
%
% Inputs:
%    signal - Kleene signal indicating the value of the four-valued signal if it is not inconclusive
%    isInconclusive - Boolean pointSegmentSignal indicating where the four-valued signal is inconclusive
%
% Outputs:
%    sig - generated fourValuedSignal object
%
% Example:
%    signal = kleeneSignal.indicator(stlInterval(1,3),kleene.True,kleene.Unknown);
%    signal = signal & kleeneSignal.indicator(stlInterval(4,6),kleene.False,kleene.True);
%    isInconclusive = pointSegmentSignal.indicator(stlInterval(2,5),true,false);
%    sig = fourValuedSignal(signal,isInconclusive);
%    sig.plot()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: logicSignal, pointSegmentSignal, kleeneSignal

% Authors:       Florian Lercher
% Written:       14-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    signal kleeneSignal
    isInconclusive pointSegmentSignal
end

methods
    % constructor
    function obj = fourValuedSignal(signal,isInconclusive)
        obj.signal = signal;
        obj.isInconclusive = isInconclusive;
    end

    % ------basic operations------
    % check if two four-valued signals are equal
    function equal = eq(obj,other)
        if isa(other,'fourValuedSignal')
            equal = obj.signal == other.signal && obj.isInconclusive == other.isInconclusive;
        else
            equal = false;
        end
    end

    % query the value of the signal at a given time
    function val = at(obj,time)
        if obj.isInconclusive.at(time)
            val = fourValued.Inconclusive;
        else
            val = fourValued.fromKleene(obj.signal.at(time));
        end
    end

    % set the signal value in a given time interval
    function sig = set(obj,interval,value)
        [kleeneVal,inc] = fourValuedSignal.u2ToKleeneInc(value);
        newSignal = obj.signal.set(interval,kleeneVal);
        newInc = obj.isInconclusive.set(interval,inc);
        sig = fourValuedSignal(newSignal,newInc);
    end

    % convert to a Kleene signal
    % by setting mapInconclusiveTo to kleene.True, we obtain the maximal refinement,
    % by setting it to kleene.False, we obtain the minimal refinement
    function sig = toKleeneSignal(obj,mapInconclusiveTo)
        if nargin < 2
            if obj.isInconclusive.anyTrue(stlInterval(0,inf))
                throw(CORAerror('CORA:notDefined','Cannot convert partial Kleene signal with inconclusive values to Kleene signal'));
            end
            sig = obj.signal;
        else
            intervals = obj.isInconclusive.findIntervals();
            sig = obj.signal;
            for i = 1:length(intervals)
                sig = sig.set(intervals(i),mapInconclusiveTo);
            end
        end
    end

    % find the time intervals in which the signal has the value match
    function intervals = findIntervals(obj,match)
        if nargin < 2
            match = fourValued.True;
        end
        switch match
            case fourValued.Inconclusive
                intervals = obj.isInconclusive.findIntervals();
            case fourValued.True
                intervals = obj.toKleeneSignal(kleene.False).findIntervals(kleene.True);
            case fourValued.False
                intervals = obj.toKleeneSignal(kleene.True).findIntervals(kleene.False);
            case fourValued.Unknown
                intervals = obj.toKleeneSignal(kleene.True).findIntervals(kleene.Unknown);
            otherwise
                assert(false,'Not a known U2 value');
        end
    end

    % ------logical operators------
    % logical negation
    function sig = not(obj)
        sig = fourValuedSignal(~obj.signal,obj.isInconclusive);
    end

    % until combination of two signals
    function sig = until(lhs,interval,rhs)
        lhsMinRefinement = lhs.toKleeneSignal(kleene.False);
        rhsMinRefinement = rhs.toKleeneSignal(kleene.False);
        minUntil = lhsMinRefinement.until(interval,rhsMinRefinement);

        lhsMaxRefinement = lhs.toKleeneSignal(kleene.True);
        rhsMaxRefinement = rhs.toKleeneSignal(kleene.True);
        maxUntil = lhsMaxRefinement.until(interval,rhsMaxRefinement);
        
        inc = kleeneSignal.combine(@(args) kleene.fromBool(args(1) ~= args(2)),minUntil,maxUntil).toBoolSignal();

        sig = fourValuedSignal(minUntil,inc);
    end

    % ------operations for treating inconclusive as undefined------
    % merge the signal with another signal
    % obj is set to the value of other wherever other is not inconclusive
    function sig = merge(obj,other)
        sig = obj;
        
        trueInt = other.findIntervals(fourValued.True);
        for i = 1:length(trueInt)
            sig = sig.set(trueInt(i),fourValued.True);
        end

        falseInt = other.findIntervals(fourValued.False);
        for i = 1:length(falseInt)
            sig = sig.set(falseInt(i),fourValued.False);
        end

        unknownInt = other.findIntervals(fourValued.Unknown);
        for i = 1:length(unknownInt)
            sig = sig.set(unknownInt(i),fourValued.Unknown);
        end
    end

    % find the largest time interval I containing 0
    % so that the signal is conclusive throughout I
    % ignoring the values in the interval ignore
    function int = conclusiveInterval(obj,ignore)
        isConclusive = ~(obj.isInconclusive.set(ignore,false));
        if ~isConclusive.at(0)
            int = stlInterval();
        else
            [time,inclusive] = isConclusive.cutAtFirstFallingEdge().findLast(true);
            int = stlInterval(0,time,true,inclusive);
        end
    end

    % ------plotting------
    % plot the signal
    function han = plot(obj)
        plotSig = obj.plotSig();
        if nargout > 0
            han = plotSig.plot();
        else
            plotSig.plot();
        end
    end

    % convert the signal to a double signal for plotting
    function sig = plotSig(obj)
        incDouble = pointSegmentSignal.combine(@double,obj.isInconclusive);
        sig = pointSegmentSignal.combine(...
            @(args) double(fourValuedSignal.kleeneIncToU2(kleene(args(1)),logical(args(2)))),...
            obj.signal.plotSig(),incDouble);
    end
end

methods (Static)
    % ------constructor-like methods------
    % indicator signals
    function sig = indicator(interval,val,default)
        [valKleene,valUnk] = fourValuedSignal.u2ToKleeneInc(val);
        [defaultKleene,defaultUnk] = fourValuedSignal.u2ToKleeneInc(default);
        kleeneSig = kleeneSignal.indicator(interval,valKleene,defaultKleene);
        unknown = pointSegmentSignal.indicator(interval,valUnk,defaultUnk);
        sig = fourValuedSignal(kleeneSig,unknown);
    end

    % embedding of Kleene signals
    function sig = fromKleeneSignal(kSig)
        inc = pointSegmentSignal.indicator(stlInterval(),false,false);
        sig = fourValuedSignal(kSig,inc);
    end

    % create a signal with a constant value
    function sig = uniformSignal(val)
        sig = fourValuedSignal.indicator(stlInterval(),val,val);
    end

    % ------point-wise operations------
    % arbitrary operation specified by op working on U2 values
    function sigs = combine(op,varargin)
        arr = [varargin{:}];
        inconclusiveKleene = arrayfun(@(s) kleeneSignal.fromBoolSignal(s.isInconclusive),arr);
        c = num2cell([arr.signal,inconclusiveKleene]);
        pSigs = kleeneSignal.combine(@u2Op,c{:});
        for i = 1:2:size(pSigs,2)
            sigs(ceil(i / 2)) = fourValuedSignal(pSigs(i),pSigs(i + 1).toBoolSignal());
        end

        % adapt op to work with our representation
        function out = u2Op(args)
            u2s = arrayfun(@fourValuedSignal.kleeneIncToU2,args(1:end / 2),args(end / 2 + 1:end));
            res = op(u2s);
            for j = 1:2:length(res) * 2
                [boolVal,inc] = fourValuedSignal.u2ToKleeneInc(res(ceil(j / 2)));
                out(j) = boolVal;
                out(j + 1) = kleene.fromBool(inc);
            end
        end
    end

    % logical conjunction
    function sig = and_(varargin)
        minRefinements = cellfun(@(s) s.toKleeneSignal(kleene.False),varargin,'UniformOutput',false);
        minAnd = kleeneSignal.and_(minRefinements{:});

        maxRefinements = cellfun(@(s) s.toKleeneSignal(kleene.True),varargin,'UniformOutput',false);
        maxAnd = kleeneSignal.and_(maxRefinements{:});

        inc = kleeneSignal.combine(@(args) kleene.fromBool(range(args) ~= 0),minAnd,maxAnd).toBoolSignal();

        sig = fourValuedSignal(minAnd,inc);
    end

    % logical disjunction
    function sig = or_(varargin)
        minRefinements = cellfun(@(s) s.toKleeneSignal(kleene.False),varargin,'UniformOutput',false);
        minOr = kleeneSignal.or_(minRefinements{:});

        maxRefinements = cellfun(@(s) s.toKleeneSignal(kleene.True),varargin,'UniformOutput',false);
        maxOr = kleeneSignal.or_(maxRefinements{:});

        inc = kleeneSignal.combine(@(args) kleene.fromBool(range(args) ~= 0),minOr,maxOr).toBoolSignal();

        sig = fourValuedSignal(minOr,inc);
    end
end

methods (Static,Access = private)
    % ------helper functions for representation------
    % from representation to U2 value
    function u = kleeneIncToU2(kleeneVal,inc)
        if islogical(inc) && inc || isa(inc,'kleene') && inc == kleene.True
            u = fourValued.Inconclusive;
        else
            u = fourValued.fromKleene(kleeneVal);
        end
    end

    % from U2 value to representation
    function [kleeneVal,inc] = u2ToKleeneInc(u)
        switch u
            case fourValued.True
                kleeneVal = kleene.True;
                inc = false;
            case fourValued.False
                kleeneVal = kleene.False;
                inc = false;
            case fourValued.Unknown
                kleeneVal = kleene.Unknown;
                inc = false;
            case fourValued.Inconclusive
                kleeneVal = kleene.Unknown; % exact value does not matter
                inc = true;
            otherwise
                assert(false,'Not a known U2 value');
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
