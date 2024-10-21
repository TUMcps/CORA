classdef kleeneSignal < logicSignal
% kleeneSignal - a signal over Kleene values represented by two Boolean signals
%
% Since the Kleene signal is represented by pointSegmentSignal objects,
% most operations are simple wrappers around those of pointSegmentSignal.
%
% Syntax:
%    sig = kleeneSignal(signal,isUnknown)
%
% Inputs:
%    signal - Boolean pointSegmentSignal indicating the value of the Kleene signal if it is not unknown
%    isUnknown - Boolean pointSegmentSignal indicating where the Kleene signal is unknown
%
% Outputs:
%    sig - generated kleeneSignal object
%
% Example:
%    signal = pointSegmentSignal.indicator(stlInterval(1,3),true,false);
%    isUnknown = pointSegmentSignal.indicator(stlInterval(2,4),true,false);
%    sig = kleeneSignal(signal,isUnknown);
%    sig.plot()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: logicSignal, pointSegmentSignal

% Authors:       Florian Lercher
% Written:       14-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    signal pointSegmentSignal
    isUnknown pointSegmentSignal
end

methods
    % constructor
    function obj = kleeneSignal(signal,isUnknown)
        obj.signal = signal;
        obj.isUnknown = isUnknown;
    end

    % ------basic operations------
    % check if two Kleene signals are equal
    function equal = eq(obj,other)
        if isa(other,'kleeneSignal')
            equal = obj.signal == other.signal && obj.isUnknown == other.isUnknown;
        else
            equal = false;
        end
    end
    
    % query the value of the signal at a given time
    function val = at(obj,time)
        if obj.isUnknown.at(time)
            val = kleene.Unknown;
        else
            val = kleene.fromBool(obj.signal.at(time));
        end
    end

    % set the signal value in a given time interval
    function sig = set(obj,interval,value)
        [boolVal,unk] = kleeneSignal.kleeneToBoolUnk(value);
        val = obj.signal.set(interval,boolVal);
        unknown = obj.isUnknown.set(interval,unk);
        sig = kleeneSignal(val,unknown);
    end

    % convert to a Boolean signal
    % by setting mapUnknownTo to true, we obtain the maximal refinement,
    % by setting it to false, we obtain the minimal refinement
    function sig = toBoolSignal(obj,mapUnknownTo)
        if nargin < 2
            if obj.isUnknown.anyTrue(stlInterval(0,inf))
                throw(CORAerror('CORA:notDefined','Cannot convert Kleene signal with unknown values to Boolean signal'));
            end
            sig = obj.signal;
        else
            intervals = obj.isUnknown.findIntervals();
            sig = obj.signal;
            for i = 1:length(intervals)
                sig = sig.set(intervals(i),mapUnknownTo);
            end
        end
    end

    % find the time intervals in which the signal has the value match
    function intervals = findIntervals(obj,match)
        if nargin < 2
            match = kleene.True;
        end
        switch match
            case kleene.Unknown
                intervals = obj.isUnknown.findIntervals();
            case kleene.True
                intervals = obj.toBoolSignal(false).findIntervals();
            case kleene.False
                intervals = obj.toBoolSignal(true).findIntervals(@(x) ~x);
            otherwise
                assert(false,'Not a known Kleene value');
        end
    end

    % ------logical operators------
    % logical negation
    function sig = not(obj)
        sig = kleeneSignal(~obj.signal,obj.isUnknown);
    end

    % until combination of two signals
    function sig = until(lhs,interval,rhs)
        lhsMinRefinement = lhs.toBoolSignal(false);
        rhsMinRefinement = rhs.toBoolSignal(false);
        minUntil = lhsMinRefinement.until(interval,rhsMinRefinement);

        lhsMaxRefinement = lhs.toBoolSignal(true);
        rhsMaxRefinement = rhs.toBoolSignal(true);
        maxUntil = lhsMaxRefinement.until(interval,rhsMaxRefinement);

        unk = pointSegmentSignal.combine(@(args) args(1) ~= args(2),minUntil,maxUntil);

        sig = kleeneSignal(minUntil,unk);
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
        sig = pointSegmentSignal.combine(...
            @(args) double(kleeneSignal.boolUnkToKleene(args(1),args(2))),...
            obj.signal,obj.isUnknown);
    end
end

methods (Static)
    % ------constructor-like methods------
    % indicator signals
    function sig = indicator(interval,val,default)
        [valBool,valUnk] = kleeneSignal.kleeneToBoolUnk(val);
        [defaultBool,defaultUnk] = kleeneSignal.kleeneToBoolUnk(default);
        boolSig = pointSegmentSignal.indicator(interval,valBool,defaultBool);
        unknown = pointSegmentSignal.indicator(interval,valUnk,defaultUnk);
        sig = kleeneSignal(boolSig,unknown);
    end

    % embedding of Boolean signals
    function sig = fromBoolSignal(bSig)
        unk = pointSegmentSignal.indicator(stlInterval(),false,false);
        sig = kleeneSignal(bSig,unk);
    end

    % ------point-wise operations------
    % arbitrary operation specified by op working on Kleene values
    function sigs = combine(op,varargin)
        arr = [varargin{:}];
        c = num2cell([arr.signal,arr.isUnknown]);
        pSigs = pointSegmentSignal.combine(@kOp,c{:});
        for i = 1:2:size(pSigs,2)
            sigs(ceil(i / 2)) = kleeneSignal(pSigs(i),pSigs(i + 1));
        end

        % adapt op to work with our representation
        function out = kOp(args)
            kleenes = arrayfun(@kleeneSignal.boolUnkToKleene,args(1:end / 2),args(end / 2 + 1:end));
            res = op(kleenes);
            for j = 1:2:length(res) * 2
                [boolVal,unk] = kleeneSignal.kleeneToBoolUnk(res(ceil(j / 2)));
                out(j) = boolVal;
                out(j + 1) = unk;
            end
        end
    end

    % logical conjunction
    function sig = and_(varargin)
        minRefinements = cellfun(@(s) s.toBoolSignal(false),varargin,'UniformOutput',false);
        minAnd = pointSegmentSignal.and_(minRefinements{:});
        
        maxRefinements = cellfun(@(s) s.toBoolSignal(true),varargin,'UniformOutput',false);
        maxAnd = pointSegmentSignal.and_(maxRefinements{:});

        unk = pointSegmentSignal.combine(@(args) range(args) ~= 0,minAnd,maxAnd);

        sig = kleeneSignal(minAnd,unk);
    end

    % logical disjunction
    function sig = or_(varargin)
        minRefinements = cellfun(@(s) s.toBoolSignal(false),varargin,'UniformOutput',false);
        minOr = pointSegmentSignal.or_(minRefinements{:});
        
        maxRefinements = cellfun(@(s) s.toBoolSignal(true),varargin,'UniformOutput',false);
        maxOr = pointSegmentSignal.or_(maxRefinements{:});

        unk = pointSegmentSignal.combine(@(args) range(args) ~= 0,minOr,maxOr);

        sig = kleeneSignal(minOr,unk);
    end
end

methods (Static,Access = private)
    % ------helper functions for representation------
    % from representation to Kleene value
    function k = boolUnkToKleene(bool,unk)
        if unk
            k = kleene.Unknown;
        else
            k = kleene.fromBool(bool);
        end
    end

    % from Kleene value to representation
    function [boolVal,unk] = kleeneToBoolUnk(k)
        switch k
            case kleene.True
                boolVal = true;
                unk = false;
            case kleene.False
                boolVal = false;
                unk = false;
            case kleene.Unknown
                boolVal = false; % exact value does not matter
                unk = true;
            otherwise
                assert(false,'Not a known Kleene value');
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
