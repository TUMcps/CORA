classdef simResultAnalyzer < handle
% simResultAnalyzer - construct boolean signals for the given atomic
% propositions over simulation results
%
% Syntax:
%    analyzer = simResultAnalyzer(aps, 10)
%
% Inputs:
%    aps - array of sets representing the atomic propositions
%    dur - duration of the signal
%    index - index of the simulation to analyze
%
% Outputs:
%    obj - simResultAnalyzer object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Benedikt Seidl
% Written:       04-January-2022
% Last update:   04-July-2024 (TL, updated analyze() due to new simResult)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % map of propositions to their geometric representation
    aps containers.Map

    % the simulation result to analyze
    simRes
end

methods

    function obj = simResultAnalyzer(aps, simRes)
        obj.aps = aps;
        obj.simRes = simRes;
    end

    function sigs = analyze(obj, i, j)
        % analyze all time intervals of simulation result with index n
        k = keys(obj.aps);
        sigs = containers.Map;
        dur = obj.simRes(i).t{j}(end);

        % prepare signals for all atomic propositions
        for m = 1:length(k)
            sigs(k{m}) = signal(dur, false);
        end

        for l = 1:length(obj.simRes(i).t{j})-1
            % obtain value of simulation
            %
            % Taking the first point only is an approximation.
            % However, since simulation results have no claim of
            % correctness this inaccuracy is acceptable here.
            point = obj.simRes(i).x{j}(l,:);

            % set location property
            if ~isempty(obj.simRes(i).loc)
                loc = obj.simRes(i).loc(j);
            else
                loc = 1;
            end

            % calculate time interval
            from = obj.simRes(i).t{j}(l);
            to = obj.simRes(i).t{j}(l+1);

            % check all atomic propositions
            for m = 1:length(k)
                if obj.aps(k{m}).evaluatePoint(point.', loc)
                    sig = signal.indicator(dur, interval(from, to), true);

                    sigs(k{m}) = sigs(k{m}) | sig;
                end
            end
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
