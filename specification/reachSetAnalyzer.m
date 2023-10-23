classdef reachSetAnalyzer < matlab.mixin.Copyable
% reachSetAnalyzer - construct three-valued signals for the given atomic
% propositions over reachable sets
%
% Syntax:
%    analyzer = reachSetAnalyzer(aps, 10)
%
% Inputs:
%    aps - array of sets representing the atomic propositions
%    dur - duration of the signal
%    masks - mask signals for all atomic propositions
%
% Outputs:
%    obj - reachSetAnalyzer object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Benedikt Seidl
% Written:       24-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % map of propositions to their geometric representation
    aps containers.Map

    % Kleene signal builders for all atomic propositions
    sigBuilder containers.Map

    % mask signals for every atomic proposition
    masks containers.Map
end

methods

    function obj = reachSetAnalyzer(aps, dur, masks)
        obj.aps = aps;
        obj.sigBuilder = containers.Map;

        if nargin < 3
            obj.masks = containers.Map;
        else
            obj.masks = masks;
        end

        k = keys(aps);

        for i = 1:length(k)
            obj.sigBuilder(k{i}) = kleeneSignalBuilder(dur);
        end
    end

    function out = analyze(obj, R, i)
        % recursively analyze the list of reachable sets
        if nargin > 2
            obj.analyzeSingle(R(i))
        else
            i = 0;
        end

        % find all children of the current node
        c = children(R,i);

        if isempty(c)
            out = {obj.signals()};
        else
            % analyze the child branches
            out = {};

            for j=1:length(c)-1
                % copy this object for independent analysis
                out = [out copy(obj).analyze(R,c(j))];
            end

            % use the current copy for last child
            out = [out obj.analyze(R,c(1))];
        end
    end

    function analyzeSingle(obj, R)
        % analyze all time intervals of a single reachable set
        for i = 1:length(R.timeInterval.time)
            set = R.timeInterval.set{i};
            int = R.timeInterval.time{i};

            % we convert set into zonotope in case it was a zonoBundle
            if isa(set, 'zonoBundle')
                set = zonotope(set);
            end

            obj.analyzeTimeInterval(set, R.loc, int);
        end
    end

    function analyzeTimeInterval(obj, set, loc, int)
        % iterate through all properties
        k = keys(obj.aps);

        for i = 1:length(k)
            % get objects from map
            ap = obj.aps(k{i});
            builder = obj.sigBuilder(k{i});

            % check mask
            if ~isKey(obj.masks, k{i}) || max(obj.masks(k{i}), int)
                % can set be inside of atomic proposition?
                if ap.canBeTrue(set, loc)
                    builder.setPossiblyTrue(int);
                end

                % can set be outside of atomic proposition?
                if ap.canBeFalse(set, loc)
                    builder.setPossiblyFalse(int);
                end
            end
        end
    end

    function out = signals(obj)
        % construct three-valued signals for all properties
        out = containers.Map;

        k = keys(obj.aps);

        for i = 1:length(k)
            out(k{i}) = obj.sigBuilder(k{i}).kleeneSignal();
        end
    end

end

methods (Access = protected)

    function cp = copyElement(obj)
        % shallow copy object
        cp = copyElement@matlab.mixin.Copyable(obj);

        % create new map
        k = keys(cp.aps);
        cp.sigBuilder = containers.Map;

        % copy all signal builders
        for i = 1:length(k)
            cp.sigBuilder(k{i}) = copy(obj.sigBuilder(k{i}));
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
